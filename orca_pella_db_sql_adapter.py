import json
import psycopg2


class OrcaPellaSQLAdapter:
    _username = _password = _db_name = _host = None
    _connection = None

    def __init__(self):
        raise NotImplementedError(
            "This is strictly a static method oriented class with class variables"
        )

    @staticmethod
    def _raise_no_details_error():
        raise PermissionError(
            'use "OrcaPellaSQLAdapter.set_connection_info" first to set details'
        )

    @staticmethod
    def _get_connection_from_details(
        username: str, password: str, db_name: str, host: str = "localhost"
    ):
        conn = psycopg2.connect(
            host=host, database=db_name, user=username, password=password
        )
        return conn

    @classmethod
    def set_connection_info(
        cls, username: str, password: str, db_name: str, host: str = "localhost"
    ):
        cls._username, cls._password = username, password
        cls._db_name, cls._host = db_name, host
        cls._connection = OrcaPellaSQLAdapter._get_connection_from_details(
            username, password, db_name, host
        )

    @classmethod
    def _get_connection(cls):
        if cls._connection is None:
            OrcaPellaSQLAdapter._raise_no_details_error()
        return cls._connection

    @staticmethod
    def upload_json(filepath: str, orca_pella_schema: dict):
        conn = OrcaPellaSQLAdapter._get_connection()
        elements = orca_pella_schema["elements"]
        cursor = conn.cursor()
        json_data = json.dumps(orca_pella_schema)
        insert_query = "INSERT INTO orca_json_dump (filename, elements, json_data) VALUES (%s, %s, %s)"
        cursor.execute(insert_query, (filepath, elements, json_data))
        conn.commit()
        cursor.close()

    @classmethod
    def close_connection(cls):
        cls._connection = None

    @classmethod
    def reset_connection(cls):
        if any(
            [
                detail is None
                for detail in (cls._username, cls._password, cls._db_name, cls._host)
            ]
        ):
            OrcaPellaSQLAdapter._raise_no_details_error()
        cls._connection = OrcaPellaSQLAdapter._get_connection_from_details(
            cls._username, cls._password, cls._db_name, cls._host
        )
        return cls._connection
