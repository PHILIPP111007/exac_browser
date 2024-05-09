import os
import pymongo

from modules.settings import settings
from modules.fastapi_globals import g


def connect_db():
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=settings.DB_HOST, port=settings.DB_PORT)
    return client[settings.DB_NAME]


def get_db():
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    g.db_conn = connect_db()
    return g.db_conn
