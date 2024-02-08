class PGRouter(object):
    """
    A router to control all database operations on models in the
    pgmols application.
    """
    def db_for_read(self, model, **hints):
        """
        Attempts to read pgmols models go to pgmols_db.
        """
        if model._meta.app_label == 'pgmols':
            return 'pgmols_db'
        return None

    def db_for_write(self, model, **hints):
        """
        Attempts to write pgmols models go to pgmols_db.
        """
        if model._meta.app_label == 'pgmols':
            return 'pgmols_db'
        return None

    def allow_relation(self, obj1, obj2, **hints):
        """
        Allow relations if a model in the pgmols app is involved.
        """
        if obj1._meta.app_label == 'pgmols' or \
           obj2._meta.app_label == 'pgmols':
           return True
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """
        Make sure the pgmols app only appears in the 'pgmols_db'
        database.
        """
        if app_label == 'pgmols':
            return db == 'pgmols_db'
        return None