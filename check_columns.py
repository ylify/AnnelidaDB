import os
import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'db_project.settings')
django.setup()
from django.template.loader import render_to_string
from api.models import Genome
from html.parser import HTMLParser

html = render_to_string('table_page.html', {
    'page_title': 'Genome',
    'data_type': 'genome',
    'data': Genome.objects.select_related('taxonomy').all()[:1],
    'page_summary_cards': []
})

class TableCounter(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_thead = False
        self.in_tbody = False
        self.in_tr = False
        self.header_count = 0
        self.row_count = 0
        self.row_done = False

    def handle_starttag(self, tag, attrs):
        if tag == 'thead':
            self.in_thead = True
        elif tag == 'tbody':
            self.in_tbody = True
        elif tag == 'tr' and self.in_thead:
            self.in_tr = True
        elif tag == 'tr' and self.in_tbody and not self.row_done:
            self.in_tr = True
        elif tag in ('td', 'th') and self.in_tr:
            if self.in_thead:
                self.header_count += 1
            elif self.in_tbody:
                self.row_count += 1

    def handle_endtag(self, tag):
        if tag == 'thead':
            self.in_thead = False
        elif tag == 'tbody':
            self.in_tbody = False
        elif tag == 'tr' and self.in_tr:
            if self.in_tbody and not self.row_done:
                self.row_done = True
            self.in_tr = False

parser = TableCounter()
parser.feed(html)
print('header_count', parser.header_count)
print('row_count', parser.row_count)
print(html.split('<tbody')[1][:400])
