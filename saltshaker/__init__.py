"""SaltShaker: Pattern classification and visualization for MitoSAlt"""

from .config import ClassificationConfig
from .event_caller import EventCaller
from .spatial import SpatialGroupAnalyzer

__version__ = "0.1.0"
__all__ = ["ClassificationConfig", "EventCaller", "SpatialGroupAnalyzer"]