import json
import urllib.error
import urllib.parse
import urllib.request
from typing import Iterable, List, Optional

from django.conf import settings
from django.core.mail.backends.base import BaseEmailBackend
from django.core.mail.backends.smtp import EmailBackend as SMTPEmailBackend
from django.core.mail.message import EmailMessage

try:
    import msal
except ImportError:
    msal = None

_GRAPH_MSAL_APP = None


def is_graph_email_enabled() -> bool:
    return bool(
        msal
        and getattr(settings, 'GRAPH_CLIENT_ID', None)
        and getattr(settings, 'GRAPH_TENANT_ID', None)
        and getattr(settings, 'GRAPH_CLIENT_SECRET', None)
        and getattr(settings, 'GRAPH_MAIL_SENDER', None)
    )


def _get_msal_app():
    global _GRAPH_MSAL_APP
    if _GRAPH_MSAL_APP is None:
        if not msal:
            raise RuntimeError('msal is required for Graph mail delivery')
        authority = f"https://login.microsoftonline.com/{settings.GRAPH_TENANT_ID}"
        _GRAPH_MSAL_APP = msal.ConfidentialClientApplication(
            client_id=settings.GRAPH_CLIENT_ID,
            authority=authority,
            client_credential=settings.GRAPH_CLIENT_SECRET,
        )
    return _GRAPH_MSAL_APP


def _get_graph_access_token() -> Optional[str]:
    if not is_graph_email_enabled():
        return None

    app = _get_msal_app()
    result = app.acquire_token_silent(['https://graph.microsoft.com/.default'], account=None)
    if not result or 'access_token' not in result:
        result = app.acquire_token_for_client(scopes=['https://graph.microsoft.com/.default'])

    if not result:
        return None
    if 'access_token' in result:
        return result['access_token']

    error = result.get('error')
    description = result.get('error_description', '')
    print(f'Microsoft Graph token error: {error} {description}')
    return None


def _format_recipients(recipients: Iterable[str]) -> List[dict]:
    return [{'emailAddress': {'address': recipient}} for recipient in recipients if recipient]


def send_graph_email(
    subject: str,
    body: str,
    to_recipients: Iterable[str],
    content_type: str = 'Text',
    cc_recipients: Optional[Iterable[str]] = None,
    bcc_recipients: Optional[Iterable[str]] = None,
) -> bool:
    token = _get_graph_access_token()
    if not token:
        raise RuntimeError('Unable to acquire Microsoft Graph access token')

    if isinstance(to_recipients, str):
        to_recipients = [to_recipients]

    payload = {
        'message': {
            'subject': subject,
            'body': {
                'contentType': content_type,
                'content': body,
            },
            'toRecipients': _format_recipients(to_recipients),
            'saveToSentItems': True,
        }
    }

    if cc_recipients:
        payload['message']['ccRecipients'] = _format_recipients(cc_recipients)
    if bcc_recipients:
        payload['message']['bccRecipients'] = _format_recipients(bcc_recipients)

    url = f"https://graph.microsoft.com/v1.0/users/{urllib.parse.quote(settings.GRAPH_MAIL_SENDER)}/sendMail"
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode('utf-8'),
        headers={
            'Authorization': f'Bearer {token}',
            'Content-Type': 'application/json',
        },
        method='POST',
    )

    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            return response.status in {200, 202}
    except urllib.error.HTTPError as exc:
        error_body = exc.read().decode('utf-8', errors='ignore')
        print(f'Microsoft Graph sendMail HTTP error: {exc.code} {exc.reason} {error_body}')
        raise
    except urllib.error.URLError as exc:
        print(f'Microsoft Graph sendMail URL error: {exc}')
        raise


class GraphEmailBackend(BaseEmailBackend):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.smtp_backend = SMTPEmailBackend(*args, **kwargs)

    def send_messages(self, email_messages):
        if not email_messages:
            return 0

        if not is_graph_email_enabled():
            return self.smtp_backend.send_messages(email_messages)

        sent_count = 0
        for message in email_messages:
            try:
                all_to = message.to or []
                if not all_to:
                    continue

                content_type = 'Html' if getattr(message, 'content_subtype', 'plain') == 'html' else 'Text'
                body = message.body or ''
                cc = getattr(message, 'cc', []) or []
                bcc = getattr(message, 'bcc', []) or []

                send_graph_email(
                    subject=message.subject or '',
                    body=body,
                    to_recipients=all_to,
                    cc_recipients=cc,
                    bcc_recipients=bcc,
                    content_type=content_type,
                )
                sent_count += 1
            except Exception:
                if not self.fail_silently:
                    raise
        return sent_count
