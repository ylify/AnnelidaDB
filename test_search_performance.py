#!/usr/bin/env python3
"""
Gene Search Performance Test Script
Tests the performance of gene search functionality
"""
import os
import sys
import django
import time
from django.conf import settings

# Setup Django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'db_project.settings')
django.setup()

from api.models import EggnogAnnotation
from django.db import connection

def test_search_performance():
    """Test search performance with different queries"""
    test_queries = [
        'ribosomal',
        'kinase',
        'protein',
        'transferase',
        'helicase'
    ]

    print("Gene Search Performance Test")
    print("=" * 40)
    print(f"Total annotations: {EggnogAnnotation.objects.count()}")
    print()

    for query in test_queries:
        # Test icontains search
        start_time = time.time()
        results_icontains = EggnogAnnotation.objects.filter(
            description__icontains=query
        ).select_related('taxonomy')[:50]
        icontains_time = time.time() - start_time

        # Test with EXPLAIN to see if indexes are used (SQLite syntax)
        with connection.cursor() as cursor:
            cursor.execute(f"EXPLAIN QUERY PLAN SELECT * FROM api_eggnogannotation WHERE description LIKE '%{query}%' LIMIT 50")
            explain_result = cursor.fetchall()

        print(f"Query: '{query}'")
        print(f"  icontains search: {len(results_icontains)} results in {icontains_time:.4f}s")
        print("  EXPLAIN output:")
        for row in explain_result:
            print(f"    {row[0]}")
        print()

def test_index_usage():
    """Test if indexes are being used"""
    print("Index Usage Test")
    print("=" * 20)

    with connection.cursor() as cursor:
        # Check if indexes exist (SQLite syntax)
        cursor.execute("""
            SELECT name, sql
            FROM sqlite_master
            WHERE type='index' AND tbl_name='api_eggnogannotation'
            AND (name LIKE '%description%' OR name LIKE '%preferred_name%')
        """)
        indexes = cursor.fetchall()

        print("Description and preferred_name indexes:")
        for index in indexes:
            print(f"  {index[0]}: {index[1]}")

        # Check table size
        cursor.execute("SELECT COUNT(*) FROM api_eggnogannotation")
        count = cursor.fetchone()
        print(f"Table rows: {count[0]}")

        # Check if indexes are actually used
        print("\nTesting index usage with different queries:")
        test_queries = ['ribosomal', 'kinase']
        for query in test_queries:
            cursor.execute(f"EXPLAIN QUERY PLAN SELECT * FROM api_eggnogannotation WHERE description LIKE '%{query}%' LIMIT 5")
            plan = cursor.fetchall()
            print(f"  Query '{query}': {plan}")

if __name__ == '__main__':
    test_search_performance()
    test_index_usage()