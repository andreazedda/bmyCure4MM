"""
Comprehensive tests for clinic patient CRUD functionality.
Tests patient creation, editing, validation, and form errors.
"""
from __future__ import annotations

from datetime import date, timedelta
from django.contrib.auth import get_user_model
from django.test import TestCase, Client
from django.urls import reverse

from clinic.models import Patient, Assessment
from clinic.forms import PatientForm, AssessmentForm


class PatientFormValidationTests(TestCase):
    """Test patient form validations and business logic."""
    
    def test_valid_patient_form(self) -> None:
        """Test form with valid data."""
        form_data = {
            "mrn": "MM-2024-001",
            "first_name": "John",
            "last_name": "Doe",
            "birth_date": date(1960, 1, 1),
            "sex": "M",
            "diagnosis_date": date(2020, 1, 1),
            "notes": "Test patient"
        }
        form = PatientForm(data=form_data)
        self.assertTrue(form.is_valid(), f"Form errors: {form.errors}")
    
    def test_mrn_converted_to_uppercase(self) -> None:
        """Test that MRN is automatically converted to uppercase."""
        form_data = {
            "mrn": "mm-2024-001",  # lowercase
            "first_name": "John",
            "last_name": "Doe",
            "birth_date": date(1960, 1, 1),
            "sex": "M",
            "diagnosis_date": date(2020, 1, 1),
        }
        form = PatientForm(data=form_data)
        self.assertTrue(form.is_valid())
        self.assertEqual(form.cleaned_data["mrn"], "MM-2024-001")
    
    def test_duplicate_mrn_rejected(self) -> None:
        """Test that duplicate MRN is rejected."""
        # Create existing patient
        Patient.objects.create(
            mrn="MM-2024-001",
            first_name="Jane",
            last_name="Smith",
            birth_date=date(1965, 1, 1),
            sex="F",
            diagnosis_date=date(2021, 1, 1)
        )
        
        # Try to create another with same MRN
        form_data = {
            "mrn": "MM-2024-001",
            "first_name": "John",
            "last_name": "Doe",
            "birth_date": date(1960, 1, 1),
            "sex": "M",
            "diagnosis_date": date(2020, 1, 1),
        }
        form = PatientForm(data=form_data)
        self.assertFalse(form.is_valid())
        self.assertIn("mrn", form.errors)
        self.assertIn("already exists", str(form.errors["mrn"]))
    
    def test_birth_date_in_future_rejected(self) -> None:
        """Test that future birth date is rejected."""
        future_date = date.today() + timedelta(days=1)
        form_data = {
            "mrn": "MM-2024-002",
            "first_name": "John",
            "last_name": "Doe",
            "birth_date": future_date,
            "sex": "M",
            "diagnosis_date": date.today(),
        }
        form = PatientForm(data=form_data)
        self.assertFalse(form.is_valid())
        self.assertIn("birth_date", form.errors)
    
    def test_age_too_young_rejected(self) -> None:
        """Test that patient under 18 is rejected."""
        young_date = date.today() - timedelta(days=365 * 10)  # 10 years old
        form_data = {
            "mrn": "MM-2024-003",
            "first_name": "Young",
            "last_name": "Patient",
            "birth_date": young_date,
            "sex": "M",
            "diagnosis_date": date.today(),
        }
        form = PatientForm(data=form_data)
        self.assertFalse(form.is_valid())
        self.assertIn("birth_date", form.errors)
    
    def test_diagnosis_before_birth_rejected(self) -> None:
        """Test that diagnosis date before birth is rejected."""
        form_data = {
            "mrn": "MM-2024-004",
            "first_name": "John",
            "last_name": "Doe",
            "birth_date": date(1960, 1, 1),
            "sex": "M",
            "diagnosis_date": date(1959, 1, 1),  # Before birth!
        }
        form = PatientForm(data=form_data)
        self.assertFalse(form.is_valid())
        self.assertIn("Diagnosis date must be after birth date", str(form.errors))


class PatientCRUDViewTests(TestCase):
    """Test patient CRUD views."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user(
            "testuser",
            password="testpass123",
            is_staff=True
        )
        self.client = Client()
        self.client.force_login(self.user)
    
    def test_patient_list_view_loads(self) -> None:
        """Test patient list page loads."""
        response = self.client.get(reverse("clinic:patient_list"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Patient Registry")
    
    def test_patient_list_shows_welcome_banner(self) -> None:
        """Test patient list shows enhanced welcome banner."""
        response = self.client.get(reverse("clinic:patient_list"))
        self.assertContains(response, "Search & Filter")
        self.assertContains(response, "R-ISS Staging")
        self.assertContains(response, "Response Codes")
    
    def test_patient_create_form_loads(self) -> None:
        """Test patient creation form loads."""
        response = self.client.get(reverse("clinic:patient_new"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Create Patient")
    
    def test_patient_create_success(self) -> None:
        """Test successful patient creation."""
        form_data = {
            "mrn": "MM-2024-TEST",
            "first_name": "Test",
            "last_name": "Patient",
            "birth_date": "1960-01-01",
            "sex": "M",
            "diagnosis_date": "2020-01-01",
            "notes": "Test notes"
        }
        response = self.client.post(reverse("clinic:patient_new"), data=form_data)
        
        # Should redirect on success
        self.assertEqual(response.status_code, 302)
        
        # Patient should exist
        patient = Patient.objects.filter(mrn="MM-2024-TEST").first()
        self.assertIsNotNone(patient)
        self.assertEqual(patient.first_name, "Test")
    
    def test_patient_edit_view_loads(self) -> None:
        """Test patient edit form loads."""
        patient = Patient.objects.create(
            mrn="MM-EDIT-TEST",
            first_name="Edit",
            last_name="Test",
            birth_date=date(1965, 1, 1),
            sex="F",
            diagnosis_date=date(2021, 1, 1)
        )
        
        response = self.client.get(reverse("clinic:patient_edit", args=[patient.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Edit")
        self.assertContains(response, "MM-EDIT-TEST")
    
    def test_patient_detail_view_loads(self) -> None:
        """Test patient detail page loads."""
        patient = Patient.objects.create(
            mrn="MM-DETAIL-TEST",
            first_name="Detail",
            last_name="View",
            birth_date=date(1970, 1, 1),
            sex="M",
            diagnosis_date=date(2022, 1, 1)
        )
        
        response = self.client.get(reverse("clinic:patient_detail", args=[patient.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Detail View")
        self.assertContains(response, "MM-DETAIL-TEST")


class AssessmentFormValidationTests(TestCase):
    """Test assessment form validations and normal ranges."""
    
    def setUp(self) -> None:
        self.patient = Patient.objects.create(
            mrn="MM-ASSESS-TEST",
            first_name="Assess",
            last_name="Test",
            birth_date=date(1960, 1, 1),
            sex="M",
            diagnosis_date=date(2020, 1, 1)
        )
    
    def test_valid_assessment_form(self) -> None:
        """Test form with valid normal values."""
        form_data = {
            "date": date.today(),
            "m_protein_g_dl": 2.5,
            "flc_ratio": 0.8,
            "hemoglobin_g_dl": 12.5,
            "calcium_mg_dl": 9.0,
            "creatinine_mg_dl": 1.0,
            "beta2m_mg_l": 3.0,
            "albumin_g_dl": 4.0,
            "ldH_u_l": 180,
            "r_iss": "II",
            "response": "PR"
        }
        form = AssessmentForm(data=form_data)
        self.assertTrue(form.is_valid(), f"Form errors: {form.errors}")
    
    def test_assessment_date_in_future_rejected(self) -> None:
        """Test that future assessment date is rejected."""
        future_date = date.today() + timedelta(days=1)
        form_data = {
            "date": future_date,
            "m_protein_g_dl": 2.5,
        }
        form = AssessmentForm(data=form_data)
        self.assertFalse(form.is_valid())
        self.assertIn("date", form.errors)
    
    def test_very_high_m_protein_warns(self) -> None:
        """Test that very high M-Protein triggers warning."""
        form_data = {
            "date": date.today(),
            "m_protein_g_dl": 12.0,  # Very high
        }
        form = AssessmentForm(data=form_data)
        self.assertFalse(form.is_valid())
        self.assertIn("m_protein_g_dl", form.errors)
    
    def test_severe_hypercalcemia_warns(self) -> None:
        """Test that severe hypercalcemia triggers warning."""
        form_data = {
            "date": date.today(),
            "calcium_mg_dl": 15.0,  # Severe!
        }
        form = AssessmentForm(data=form_data)
        self.assertFalse(form.is_valid())
        self.assertIn("calcium_mg_dl", form.errors)
        self.assertIn("hypercalcemia", str(form.errors["calcium_mg_dl"]).lower())
    
    def test_negative_values_rejected(self) -> None:
        """Test that negative lab values are rejected."""
        form_data = {
            "date": date.today(),
            "m_protein_g_dl": -1.0,  # Negative!
        }
        form = AssessmentForm(data=form_data)
        self.assertFalse(form.is_valid())
    
    def test_normal_flc_ratio_accepted(self) -> None:
        """Test that normal FLC ratio is accepted."""
        form_data = {
            "date": date.today(),
            "flc_ratio": 0.85,  # Normal range: 0.26-1.65
        }
        form = AssessmentForm(data=form_data)
        self.assertTrue(form.is_valid())
    
    def test_extremely_abnormal_flc_warns(self) -> None:
        """Test that extremely abnormal FLC ratio triggers warning."""
        form_data = {
            "date": date.today(),
            "flc_ratio": 150.0,  # Extremely high
        }
        form = AssessmentForm(data=form_data)
        self.assertFalse(form.is_valid())
        self.assertIn("flc_ratio", form.errors)


class PatientFilteringTests(TestCase):
    """Test patient list filtering functionality."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user(
            "testuser",
            password="testpass123",
            is_staff=True
        )
        self.client = Client()
        self.client.force_login(self.user)
        
        # Create test patients
        self.patient1 = Patient.objects.create(
            mrn="MM-FILTER-1",
            first_name="John",
            last_name="Smith",
            birth_date=date(1960, 1, 1),
            sex="M",
            diagnosis_date=date(2020, 1, 1)
        )
        
        self.patient2 = Patient.objects.create(
            mrn="MM-FILTER-2",
            first_name="Jane",
            last_name="Doe",
            birth_date=date(1965, 1, 1),
            sex="F",
            diagnosis_date=date(2021, 1, 1)
        )
    
    def test_search_by_surname(self) -> None:
        """Test searching patients by surname."""
        response = self.client.get(reverse("clinic:patient_list"), {"q": "Smith"})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Smith")
        self.assertNotContains(response, "Doe")
    
    def test_search_case_insensitive(self) -> None:
        """Test search is case-insensitive."""
        response = self.client.get(reverse("clinic:patient_list"), {"q": "smith"})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Smith")


class PatientIntegrationTests(TestCase):
    """Integration tests for complete patient workflows."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user(
            "testuser",
            password="testpass123",
            is_staff=True
        )
        self.client = Client()
        self.client.force_login(self.user)
    
    def test_complete_patient_workflow(self) -> None:
        """Test complete workflow: create patient, add assessment, view detail."""
        # Step 1: Create patient
        form_data = {
            "mrn": "MM-WORKFLOW-TEST",
            "first_name": "Workflow",
            "last_name": "Test",
            "birth_date": "1960-01-01",
            "sex": "M",
            "diagnosis_date": "2020-01-01",
        }
        response = self.client.post(reverse("clinic:patient_new"), data=form_data)
        self.assertEqual(response.status_code, 302)
        
        # Step 2: Verify patient exists
        patient = Patient.objects.filter(mrn="MM-WORKFLOW-TEST").first()
        self.assertIsNotNone(patient)
        
        # Step 3: Add assessment
        assessment_data = {
            "date": date.today().isoformat(),
            "m_protein_g_dl": 3.0,
            "response": "PR"
        }
        response = self.client.post(
            reverse("clinic:assessment_new", args=[patient.pk]),
            data=assessment_data
        )
        
        # Step 4: View patient detail
        response = self.client.get(reverse("clinic:patient_detail", args=[patient.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Workflow Test")
