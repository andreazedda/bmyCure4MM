import datetime

from django.test import TestCase

from clinic import models


class PatientModelTest(TestCase):
    def test_create_patient(self):
        patient = models.Patient.objects.create(
            mrn="TEST-001",
            first_name="Test",
            last_name="Patient",
            birth_date=datetime.date(1980, 1, 1),
            sex="M",
            diagnosis_date=datetime.date(2020, 1, 1),
        )
        self.assertIsNotNone(patient.pk)
