'''Errors and exceptions for antiSMASH'''


class AntiSmashError(RuntimeError):
    '''Base class for all antiSMASH custom errors'''
    pass


class NoInputFileError(AntiSmashError):
    '''Thrown when an input file doesn't exist'''
    pass


class NCBIDownloadError(AntiSmashError):
    '''Thrown when download error messages from NCBI downloads are detected'''
    pass


class InvalidIdError(AntiSmashError):
    '''Thrown when running into invalid gene IDs during runtime'''
    pass


class InvalidLocationError(AntiSmashError):
    '''Thrown when running into invalid gene locations during runtime'''
    pass


class DuplicatePromoterError(AntiSmashError):
    '''Thrown when running into valid but duplicate promoter sequences during runtime'''
    pass
