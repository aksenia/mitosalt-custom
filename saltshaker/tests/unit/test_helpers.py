"""
Unit tests for EventCaller helper functions

Tests the static parsing methods that process comma-separated values
from cluster files.
"""
import pytest
from saltshaker.event_caller import EventCaller


class TestParseCommaSeparatedMedian:
    """Tests for _parse_comma_separated_median"""
    
    def test_single_value(self):
        """Single value should return that value as float"""
        assert EventCaller._parse_comma_separated_median("100") == 100.0
    
    def test_two_values(self):
        """Two values should return their average"""
        assert EventCaller._parse_comma_separated_median("100,200") == 150.0
    
    def test_three_values_odd(self):
        """Three values should return middle value"""
        assert EventCaller._parse_comma_separated_median("100,102,104") == 102.0
    
    def test_unsorted_values(self):
        """Values should be sorted before finding median"""
        assert EventCaller._parse_comma_separated_median("50,10,30,20,40") == 30.0
    
    def test_with_spaces(self):
        """Should handle spaces after commas"""
        result = EventCaller._parse_comma_separated_median("100, 102, 104")
        assert result == 102.0
    
    def test_large_list(self):
        """Should handle larger lists correctly"""
        values = ",".join(str(i) for i in range(1, 101))
        assert EventCaller._parse_comma_separated_median(values) == 50.5


class TestParseCommaSeparatedMin:
    """Tests for _parse_comma_separated_min"""
    
    def test_single_value(self):
        """Single value should return that value"""
        assert EventCaller._parse_comma_separated_min("100") == 100
    
    def test_multiple_values(self):
        """Should return minimum value"""
        assert EventCaller._parse_comma_separated_min("50,10,30,20,40") == 10
    
    def test_unsorted_values(self):
        """Should work with unsorted values"""
        assert EventCaller._parse_comma_separated_min("100,50,200,25") == 25


class TestParseCommaSeparatedMax:
    """Tests for _parse_comma_separated_max"""
    
    def test_single_value(self):
        """Single value should return that value"""
        assert EventCaller._parse_comma_separated_max("100") == 100
    
    def test_multiple_values(self):
        """Should return maximum value"""
        assert EventCaller._parse_comma_separated_max("50,10,30,20,40") == 50
    
    def test_unsorted_values(self):
        """Should work with unsorted values"""
        assert EventCaller._parse_comma_separated_max("100,50,200,25") == 200
