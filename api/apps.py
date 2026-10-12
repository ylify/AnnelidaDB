from django.apps import AppConfig


class ApiConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'api'

    def ready(self):
        try:
            from .views import start_statistics_refresh_thread, start_ncbi_omics_refresh_thread
            start_statistics_refresh_thread()
            start_ncbi_omics_refresh_thread()
        except Exception:
            pass
