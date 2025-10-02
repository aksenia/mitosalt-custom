"""SaltShaker: Pattern classification and visualization for MitoSAlt"""

from .config import ClassificationConfig
from .event_caller import EventCaller
from .spatial import SpatialGroupAnalyzer
from .classifier import EventClassifier
from . import utils
from .visualizer import CircularPlotter, plot_circular

__version__ = "0.1.0"
__all__ = ["ClassificationConfig", "EventCaller", "EventClassifier", "SpatialGroupAnalyzer", "utils", "CircularPlotter", "plot_circular"]