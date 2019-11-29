class EmptyFile(Exception):
    def __init__(self, message="File is empty!"):
        super().__init__(message)