from django.db import migrations


def seed_and_guard(apps, schema_editor):
    apps.get_model("accounts", "AccessPolicy").objects.get_or_create(pk=1)
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("SELECT COUNT(*) FROM (SELECT lower(trim(email)) FROM auth_user WHERE trim(email) <> '' GROUP BY lower(trim(email)) HAVING COUNT(*) > 1) conflicts")
        if cursor.fetchone()[0]:
            raise RuntimeError("Duplicate legacy email identities require resolution; no accounts were merged.")
        cursor.execute("CREATE UNIQUE INDEX bmy_auth_email_identity ON auth_user (lower(trim(email))) WHERE trim(email) <> ''")


def reverse_guard(apps, schema_editor):
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("DROP INDEX IF EXISTS bmy_auth_email_identity")


class Migration(migrations.Migration):
    dependencies = [("accounts", "0001_initial"), ("auth", "0012_alter_user_first_name_max_length")]
    operations = [migrations.RunPython(seed_and_guard, reverse_guard)]
