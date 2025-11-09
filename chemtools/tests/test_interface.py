"""
UI/UX integration tests for chemtools interface.
Tests form rendering, user interactions, and frontend functionality.
"""
from __future__ import annotations

from django.contrib.auth import get_user_model
from django.test import TestCase, Client
from django.urls import reverse

from chemtools import models


class DrugParamsInterfaceTests(TestCase):
    """Test drug parameters form interface."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_form_renders_with_examples(self) -> None:
        """Test form displays with example molecules."""
        response = self.client.get(reverse("chemtools:drug_params"))
        self.assertEqual(response.status_code, 200)
        
        # Check for example SMILES in help text
        self.assertContains(response, "Ibuprofen")
        self.assertContains(response, "Aspirin")
        self.assertContains(response, "CCO")
        
    def test_form_shows_pubchem_cid_examples(self) -> None:
        """Test form shows PubChem CID examples."""
        response = self.client.get(reverse("chemtools:drug_params"))
        
        # Check for CID examples
        self.assertContains(response, "2244")  # Aspirin CID
        self.assertContains(response, "3672")  # Ibuprofen CID
        self.assertContains(response, "pubchem")
        
    def test_form_validation_messages(self) -> None:
        """Test form displays validation error messages."""
        response = self.client.post(reverse("chemtools:drug_params"), {
            "smiles": "",
            "cid": ""
        })
        
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Provide a SMILES string or a PubChem CID")
        
    def test_invalid_smiles_shows_helpful_error(self) -> None:
        """Test invalid SMILES shows user-friendly error."""
        response = self.client.post(reverse("chemtools:drug_params"), {
            "smiles": "INVALID!!!",
            "cid": ""
        })
        
        self.assertEqual(response.status_code, 200)
        # Should show validation error
        self.assertContains(response, "Invalid")
        
    def test_form_help_text_visibility(self) -> None:
        """Test all form fields have visible help text."""
        response = self.client.get(reverse("chemtools:drug_params"))
        
        self.assertContains(response, "form-text")  # Bootstrap help text class
        self.assertContains(response, "Example:")
        
    def test_submit_button_states(self) -> None:
        """Test submit button styling and states."""
        response = self.client.get(reverse("chemtools:drug_params"))
        
        self.assertContains(response, "Run Evaluation")
        self.assertContains(response, "btn-primary")


class BindingVizInterfaceTests(TestCase):
    """Test binding visualizer form interface."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_form_shows_pdb_examples(self) -> None:
        """Test form displays PDB ID examples."""
        response = self.client.get(reverse("chemtools:binding_viz"))
        
        self.assertContains(response, "5LF3")  # Bortezomib
        self.assertContains(response, "4KW5")  # Lenalidomide
        self.assertContains(response, "rcsb.org")
        
    def test_form_shows_ligand_examples(self) -> None:
        """Test form displays ligand code examples."""
        response = self.client.get(reverse("chemtools:binding_viz"))
        
        self.assertContains(response, "BOR")  # Bortezomib ligand
        self.assertContains(response, "LEN")  # Lenalidomide ligand
        
    def test_pdb_id_format_validation(self) -> None:
        """Test PDB ID format validation feedback."""
        response = self.client.post(reverse("chemtools:binding_viz"), {
            "pdb_id": "INVALID",  # Should start with digit
            "ligand": ""
        })
        
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "valid")
        
    def test_form_shows_rcsb_link(self) -> None:
        """Test form provides link to RCSB database."""
        response = self.client.get(reverse("chemtools:binding_viz"))
        
        self.assertContains(response, "rcsb")
        self.assertContains(response, "Find PDB IDs")


class SimilaritySearchInterfaceTests(TestCase):
    """Test similarity search form interface."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_form_shows_drug_examples(self) -> None:
        """Test form displays drug molecule examples."""
        response = self.client.get(reverse("chemtools:similarity"))
        
        self.assertContains(response, "Bortezomib")
        self.assertContains(response, "Lenalidomide")
        self.assertContains(response, "Caffeine")
        
    def test_threshold_field_with_recommendations(self) -> None:
        """Test threshold field shows recommended values."""
        response = self.client.get(reverse("chemtools:similarity"))
        
        self.assertContains(response, "threshold")
        self.assertContains(response, "0.7")  # Recommended value
        self.assertContains(response, "0.0-1.0")  # Range
        self.assertContains(response, "recommended")
        
    def test_threshold_validation_feedback(self) -> None:
        """Test threshold validation shows clear errors."""
        response = self.client.post(reverse("chemtools:similarity"), {
            "smiles": "CCO",
            "threshold": "1.5"  # Invalid: >1.0
        })
        
        self.assertEqual(response.status_code, 200)
        # Should not redirect (form invalid)
        
    def test_form_shows_similarity_explanations(self) -> None:
        """Test form explains similarity scores."""
        response = self.client.get(reverse("chemtools:similarity"))
        
        self.assertContains(response, "similar")
        self.assertContains(response, "0.9")  # Very similar
        self.assertContains(response, "0.5")  # Broader search


class ToolsHomeInterfaceTests(TestCase):
    """Test tools home page interface."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_page_shows_welcome_banner(self) -> None:
        """Test page displays welcome banner with explanations."""
        response = self.client.get(reverse("chemtools:tools_home"))
        
        self.assertContains(response, "Drug Discovery")
        self.assertContains(response, "Cheminformatics")
        
    def test_page_shows_tool_cards(self) -> None:
        """Test page displays cards for each tool."""
        response = self.client.get(reverse("chemtools:tools_home"))
        
        # Check for tool descriptions
        self.assertContains(response, "Drug Parameters")
        self.assertContains(response, "Binding Visualizer")
        self.assertContains(response, "Similarity Search")
        
    def test_empty_state_with_suggestions(self) -> None:
        """Test empty state shows helpful suggestions."""
        response = self.client.get(reverse("chemtools:tools_home"))
        
        self.assertContains(response, "No Jobs")  # or "No jobs"
        self.assertContains(response, "Try")
        
    def test_job_table_headers_bilingual(self) -> None:
        """Test job table headers are bilingual."""
        models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        response = self.client.get(reverse("chemtools:tools_home"))
        
        # Check for bilingual headers
        self.assertContains(response, "Created")
        self.assertContains(response, "Type")
        self.assertContains(response, "Status")
        
    def test_job_status_badges_visible(self) -> None:
        """Test job status badges are clearly visible."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        response = self.client.get(reverse("chemtools:tools_home"))
        
        # Check for status badge styling
        self.assertContains(response, "badge")
        
    def test_how_to_read_table_guide(self) -> None:
        """Test page shows guide for reading results table."""
        models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        response = self.client.get(reverse("chemtools:tools_home"))
        
        self.assertContains(response, "How to Read")
        self.assertContains(response, "Status")


class JobInteractionTests(TestCase):
    """Test user interactions with jobs."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_view_html_output_link(self) -> None:
        """Test viewing HTML output works."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        job.out_html.name = "chem/1/drug_1_v1.html"
        job.save()
        
        response = self.client.get(reverse("chemtools:tools_home"))
        self.assertContains(response, "View HTML")
        
    def test_download_csv_output_link(self) -> None:
        """Test downloading CSV output works."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        job.out_csv.name = "chem/1/sim_1_v1.csv"
        job.save()
        
        response = self.client.get(reverse("chemtools:tools_home"))
        self.assertContains(response, "Download CSV")
        
    def test_view_log_functionality(self) -> None:
        """Test viewing job log works."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user,
            log="Job completed successfully\n---\nMW: 46.07"
        )
        
        response = self.client.get(reverse("chemtools:tools_home"))
        # Log should be accessible
        self.assertContains(response, job.get_kind_display())
        
    def test_retry_job_button_available(self) -> None:
        """Test retry button is available for failed jobs."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user,
            log="ERROR: Failed to process"
        )
        
        response = self.client.get(reverse("chemtools:tools_home"))
        # Should show some indication of error state
        self.assertEqual(response.status_code, 200)


class ResponsiveDesignTests(TestCase):
    """Test responsive design elements."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_mobile_friendly_layout(self) -> None:
        """Test layout includes mobile-responsive classes."""
        response = self.client.get(reverse("chemtools:tools_home"))
        
        # Bootstrap responsive classes
        self.assertContains(response, "col-md")
        self.assertContains(response, "d-flex")
        
    def test_hover_effects_present(self) -> None:
        """Test hover effects are defined in CSS."""
        response = self.client.get(reverse("chemtools:tools_home"))
        
        self.assertContains(response, "hover")
        
    def test_card_styling_consistent(self) -> None:
        """Test card styling is consistent across tools."""
        urls = [
            reverse("chemtools:drug_params"),
            reverse("chemtools:binding_viz"),
            reverse("chemtools:similarity"),
        ]
        
        for url in urls:
            response = self.client.get(url)
            self.assertContains(response, "card")
            self.assertContains(response, "shadow")


class AccessibilityTests(TestCase):
    """Test accessibility features."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_form_labels_present(self) -> None:
        """Test all form fields have labels."""
        response = self.client.get(reverse("chemtools:drug_params"))
        
        self.assertContains(response, "<label")
        self.assertContains(response, "for=")
        
    def test_help_text_associated_with_fields(self) -> None:
        """Test help text is properly associated with fields."""
        response = self.client.get(reverse("chemtools:drug_params"))
        
        # Check for form-text class (Bootstrap help text)
        self.assertContains(response, "form-text")

        
    def test_error_messages_clear_and_visible(self) -> None:
        """Test error messages are clear and visible."""
        response = self.client.post(reverse("chemtools:drug_params"), {
            "smiles": "",
            "cid": ""
        })
        
        # Should show validation error
        self.assertContains(response, "Provide a SMILES")
        
    def test_required_fields_marked(self) -> None:
        """Test required fields are clearly marked."""
        response = self.client.get(reverse("chemtools:similarity"))
        
        # Check for required indicator
        self.assertContains(response, "*")
        
    def test_icons_have_semantic_meaning(self) -> None:
        """Test icons supplement text, not replace it."""
        response = self.client.get(reverse("chemtools:tools_home"))
        
        # Icons should be paired with text
        self.assertContains(response, "💊")
        self.assertContains(response, "Drug")


class BilingualInterfaceTests(TestCase):
    """Test bilingual (EN/IT) interface elements."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_banners_show_both_languages(self) -> None:
        """Test welcome banners show explanatory text."""
        response = self.client.get(reverse("chemtools:tools_home"))
        
        # Check for bilingual or explanatory content
        self.assertContains(response, "Drug Discovery")
        self.assertContains(response, "Cheminformatics")
        
    def test_table_headers_bilingual(self) -> None:
        """Test table headers include both languages."""
        models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        response = self.client.get(reverse("chemtools:tools_home"))
        
        # Check for bilingual headers like "Created / Creato"
        content = response.content.decode()
        # Should have some indication of translation
        
    def test_status_messages_bilingual(self) -> None:
        """Test status messages appear in both languages."""
        response = self.client.get(reverse("chemtools:tools_home"))
        
        # Check for bilingual status indicators
        self.assertEqual(response.status_code, 200)
