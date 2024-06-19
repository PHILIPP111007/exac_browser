import os

from modules.settings import settings


def get_autocomplete_strings() -> list[str]:
    if os.path.exists(settings.AUTOCOMPLETE_STRINGS_DIR):
        with open(settings.AUTOCOMPLETE_STRINGS_DIR, 'r') as f:
            return list(map(lambda x: x.strip(), f.readlines()))


AUTOCOMPLETE_STRINGS = get_autocomplete_strings()
