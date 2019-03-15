"""An abstract base class for classes in DANCE."""

from abc import ABCMeta, abstractmethod


class DanceRunBase(metaclass=ABCMeta):
    """
    DanceRunBase is an abstract base class for classes in DANCE. Classes in
    DANCE inherit from DanceRunBase if they want to use a run() method. They
    then implement the run() method.

    IMPORTANT: Though it cannot be required, derived classes should call
    check_run_fatal() at the beginning of their run() methods to check if the
    class has already been run. Also, derived classes should call the
    DanceRunBase __init__ in their __init__

    Attributes:
        _run_yet: whether the class has been run yet
    """

    def __init__(self):
        self._run_yet = False

    @abstractmethod
    def run(self):
        """A method that users call when they want the class to execute"""
        pass

    def check_run_fatal(self):
        """Checks if the class has already been run; raises RuntimeError if so"""
        if self._run_yet:
            raise RuntimeError(
                f"This {type(self).__name__} has already been run")
        self._run_yet = True
