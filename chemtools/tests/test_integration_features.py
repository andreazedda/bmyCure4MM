"""
Comprehensive tests for integrated drug discovery features.
Tests the complete platform integration including:
- Enhanced Drug Parameters with properties table
- Fixed PubChem similarity search API
- Integrated web views for all job types
- CSV parsing and HTML embedding
"""
from __future__ import annotations

import csv
import tempfile
from pathlib import Path
from unittest.mock import patch, Mock

from django.contrib.auth import get_user_model
from django.core.files.base import ContentFile
from django.test import TestCase
from django.urls import reverse

from chemtools import models


class IntegratedResultsViewTests(TestCase):
    """Test integrated results viewing for all job types."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_drug_params_shows_3d_and_properties(self) -> None:
        """Test Drug Parameters shows both 3D structure AND properties table."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        # Create realistic HTML output with both 3D viewer and properties
        html_content = """
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        </head>
        <body>
            <h2>Drug Parameter Evaluation for CCO</h2>
            
            <!-- 3D Structure Viewer -->
            <div id="container-01" style="width:600px;height:400px;"></div>
            <script>
                let viewer = $3Dmol.createViewer($("#container-01"));
                viewer.setBackgroundColor(0xffffff);
                viewer.addModel("CCO", "smi");
                viewer.setStyle({stick:{}});
                viewer.zoomTo();
                viewer.render();
            </script>
            
            <!-- Molecular Properties Table -->
            <h3>Molecular Properties</h3>
            <table class="parameters-table" style="width:100%">
                <tr>
                    <td><strong>Molecular Weight (MW)</strong></td>
                    <td>46.07 g/mol</td>
                    <td><span style="color:#28a745">✓</span> Optimal (&lt; 500)</td>
                </tr>
                <tr>
                    <td><strong>Partition Coefficient (LogP)</strong></td>
                    <td>-0.31</td>
                    <td><span style="color:#28a745">✓</span> Optimal (0-5)</td>
                </tr>
                <tr>
                    <td><strong>H-Bond Donors (HBD)</strong></td>
                    <td>1</td>
                    <td><span style="color:#28a745">✓</span> Optimal (≤ 5)</td>
                </tr>
                <tr>
                    <td><strong>H-Bond Acceptors (HBA)</strong></td>
                    <td>1</td>
                    <td><span style="color:#28a745">✓</span> Optimal (≤ 10)</td>
                </tr>
                <tr>
                    <td><strong>Topological Polar Surface Area (TPSA)</strong></td>
                    <td>20.23 Ų</td>
                    <td><span style="color:#28a745">✓</span> Optimal (&lt; 140)</td>
                </tr>
                <tr>
                    <td><strong>Rotatable Bonds</strong></td>
                    <td>0</td>
                    <td><span style="color:#28a745">✓</span> Optimal (&lt; 10)</td>
                </tr>
                <tr>
                    <td><strong>Aqueous Solubility (LogS)</strong></td>
                    <td>-0.77</td>
                    <td><span style="color:#28a745">✓</span> Optimal (&gt; -4)</td>
                </tr>
            </table>
            
            <!-- Lipinski's Rule of Five -->
            <div class="alert alert-info">
                <h4>Lipinski's Rule of Five</h4>
                <p>This compound satisfies all Lipinski criteria for oral drug-likeness.</p>
            </div>
        </body>
        </html>
        """
        
        job.out_html.save("drug_params.html", ContentFile(html_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Verify 3D viewer is present
        self.assertContains(response, "$3Dmol")
        self.assertContains(response, "container-01")
        
        # Verify properties table is present
        self.assertContains(response, "parameters-table")
        self.assertContains(response, "Molecular Weight")
        self.assertContains(response, "46.07")
        self.assertContains(response, "LogP")
        self.assertContains(response, "-0.31")
        self.assertContains(response, "H-Bond Donors")
        self.assertContains(response, "TPSA")
        self.assertContains(response, "20.23")
        
        # Verify drug-likeness indicators
        self.assertContains(response, "✓")
        self.assertContains(response, "Optimal")
        self.assertContains(response, "Lipinski")
        
    def test_similarity_search_csv_to_html_table(self) -> None:
        """Test Similarity Search CSV is rendered as interactive table."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            input_b="0.85",
            user=self.user
        )
        
        # Create realistic CSV output
        csv_content = """CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name
702,CCO,1.0,46.07,-0.31,1,1,20.23,0,-0.77,Ethanol
1031,CCCO,0.95,60.1,0.25,1,1,20.23,1,-1.01,1-Propanol
6276,CC(C)O,0.92,60.1,0.05,1,1,20.23,0,-0.89,2-Propanol
887,CC(C)(C)O,0.88,74.12,0.35,1,1,20.23,0,-1.22,tert-Butanol
6584,CCCCO,0.87,74.12,0.88,1,1,20.23,2,-1.43,1-Butanol"""
        
        job.out_csv.save("similarity.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Verify table structure
        self.assertContains(response, "Similarity Search Results")
        self.assertContains(response, "<table")
        self.assertContains(response, "Ethanol")
        self.assertContains(response, "1-Propanol")
        self.assertContains(response, "2-Propanol")
        
        # Verify similarity values are displayed
        self.assertContains(response, "1.0")
        self.assertContains(response, "0.95")
        self.assertContains(response, "0.92")
        
        # Verify molecular properties are displayed
        self.assertContains(response, "46.07")  # MW
        self.assertContains(response, "-0.31")  # LogP
        
        # Verify PubChem links
        self.assertContains(response, "pubchem.ncbi.nlm.nih.gov/compound/702")
        self.assertContains(response, "pubchem.ncbi.nlm.nih.gov/compound/1031")
        
    def test_similarity_color_coding(self) -> None:
        """Test similarity results have color-coded progress bars."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="c1ccccc1",
            input_b="0.8",
            user=self.user
        )
        
        # CSV with various similarity levels
        csv_content = """CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name
241,c1ccccc1,1.0,78.11,2.13,0,0,0.0,0,-2.64,Benzene
7501,c1ccc(C)cc1,0.95,92.14,2.73,0,0,0.0,0,-3.15,Toluene
7237,c1ccc(O)cc1,0.93,94.11,1.46,1,1,20.23,0,-1.48,Phenol
7505,c1ccc(N)cc1,0.91,93.13,0.90,1,1,26.02,0,-0.92,Aniline
7794,c1ccc(Cl)cc1,0.88,112.56,2.84,0,0,0.0,0,-2.98,Chlorobenzene"""
        
        job.out_csv.save("sim.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Verify progress bars exist
        self.assertContains(response, "progress-bar")
        
        # Verify color coding based on similarity
        # Perfect match (1.0) = green
        self.assertContains(response, "bg-success")
        # High similarity (≥0.95) = blue
        self.assertContains(response, "bg-info")
        # Good similarity (0.90-0.95) = yellow
        self.assertContains(response, "bg-warning")
        
    def test_binding_viz_embedded_viewer(self) -> None:
        """Test Binding Visualizer shows embedded 3D protein-ligand viewer."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.BIND,
            input_a="5LF3",
            input_b="BOR",
            user=self.user
        )
        
        # Create realistic HTML with 3D viewer
        html_content = """
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
        </head>
        <body>
            <h2>Protein-Ligand Binding Visualization: 5LF3 - BOR</h2>
            <div id="viewport" style="width:100%;height:600px;"></div>
            <script>
                let viewer = $3Dmol.createViewer("#viewport");
                
                // Load protein structure
                $.get('/media/pdb_cache/5lf3.pdb', function(data) {
                    viewer.addModel(data, "pdb");
                    viewer.setStyle({cartoon:{color:'spectrum'}});
                    
                    // Highlight ligand
                    viewer.setStyle({resn:'BOR'}, {stick:{colorscheme:'greenCarbon'}});
                    
                    // Add surface around ligand
                    viewer.addSurface($3Dmol.SurfaceType.VDW, {
                        opacity:0.7,
                        color:'lightblue'
                    }, {resi:['BOR']});
                    
                    viewer.zoomTo({resn:'BOR'});
                    viewer.render();
                });
            </script>
            
            <div class="info-box">
                <h3>Binding Site Information</h3>
                <p><strong>PDB ID:</strong> 5LF3</p>
                <p><strong>Ligand:</strong> BOR (Bortezomib)</p>
                <p><strong>Resolution:</strong> 1.85 Å</p>
            </div>
        </body>
        </html>
        """
        
        job.out_html.save("binding.html", ContentFile(html_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Verify 3D viewer is present
        self.assertContains(response, "$3Dmol")
        self.assertContains(response, "viewport")
        self.assertContains(response, "addModel")
        self.assertContains(response, "setStyle")
        
        # Verify metadata
        self.assertContains(response, "5LF3")
        self.assertContains(response, "BOR")
        self.assertContains(response, "Bortezomib")
        
    def test_empty_similarity_results(self) -> None:
        """Test empty similarity results show helpful message."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCCCCCCCCCCCCCCCCCCC",
            input_b="0.99",
            user=self.user
        )
        
        # CSV with only headers
        csv_content = "CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name\n"
        job.out_csv.save("empty.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Verify table structure (should still render section even if empty)
        self.assertContains(response, "job.kind == 'SIM'") or self.assertContains(response, "Ligand Similarity Search")
        
    def test_progress_tracking_display(self) -> None:
        """Test progress tracking is displayed correctly."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            input_b="0.8",
            user=self.user,
            progress_percent=65,
            progress_message="Processing compound 65 of 100"
        )
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Verify progress bar
        self.assertContains(response, "progress-bar")
        self.assertContains(response, "65%")
        self.assertContains(response, "Processing compound 65 of 100")
        
    def test_auto_refresh_for_queued_jobs(self) -> None:
        """Test auto-refresh meta tag for queued/running jobs."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user,
            progress_percent=30
        )
        # No output files = still queued
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Should have auto-refresh
        self.assertContains(response, '<meta http-equiv="refresh"')
        self.assertContains(response, 'content="5"')  # 5 second refresh
        
    def test_no_refresh_for_completed_jobs(self) -> None:
        """Test no auto-refresh for completed jobs."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a="CCO",
            user=self.user
        )
        
        # Job has output = completed
        job.out_html.save("output.html", ContentFile("<html>Done</html>"), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Should NOT have auto-refresh
        response_content = response.content.decode()
        self.assertNotIn('http-equiv="refresh"', response_content)


class PubChemAPIIntegrationTests(TestCase):
    """Test correct PubChem API usage."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        
    def test_fastsimilarity_2d_endpoint(self) -> None:
        """Test that fastsimilarity_2d endpoint is used correctly."""
        with patch("LigandSimilaritySearcher.sources.lib.similarity.requests.get") as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.json.return_value = {
                "IdentifierList": {"CID": [702, 1031, 6276]}
            }
            mock_get.return_value = mock_response
            
            from LigandSimilaritySearcher.sources.lib.similarity import search_similar_compounds
            
            result = search_similar_compounds("CCO", n_results=10)
            
            # Verify correct endpoint
            called_url = mock_get.call_args[0][0]
            self.assertIn("fastsimilarity_2d", called_url)
            self.assertIn("smiles", called_url.lower())
            self.assertIn("CCO", called_url)
            
            # Verify results returned
            self.assertEqual(len(result), 3)
            
    def test_properties_api_for_smiles(self) -> None:
        """Test properties API is used to fetch SMILES."""
        with patch("LigandSimilaritySearcher.sources.lib.similarity.requests.get") as mock_get:
            # First call: similarity search
            sim_response = Mock()
            sim_response.status_code = 200
            sim_response.json.return_value = {
                "IdentifierList": {"CID": [702]}
            }
            
            # Second call: properties
            props_response = Mock()
            props_response.status_code = 200
            props_response.json.return_value = {
                "PropertyTable": {
                    "Properties": [
                        {"CID": 702, "IsomericSMILES": "CCO"}
                    ]
                }
            }
            
            mock_get.side_effect = [sim_response, props_response]
            
            from LigandSimilaritySearcher.sources.lib.similarity import search_similar_compounds
            result = search_similar_compounds("CCO", n_results=1)
            
            # Verify two API calls
            self.assertEqual(mock_get.call_count, 2)
            
            # Verify second call is properties endpoint
            second_url = mock_get.call_args_list[1][0][0]
            self.assertIn("property", second_url.lower())
            self.assertIn("702", second_url)  # CID
            self.assertIn("IsomericSMILES", second_url)
            
    def test_retry_on_api_failure(self) -> None:
        """Test retry mechanism for API failures."""
        with patch("LigandSimilaritySearcher.sources.lib.similarity.requests.get") as mock_get:
            with patch("LigandSimilaritySearcher.sources.lib.similarity.time.sleep"):  # Skip sleep
                # First call fails, second succeeds
                fail_response = Mock()
                fail_response.status_code = 500
                fail_response.raise_for_status.side_effect = Exception("Server Error")
                
                success_response = Mock()
                success_response.status_code = 200
                success_response.json.return_value = {
                    "IdentifierList": {"CID": [702]}
                }
                
                mock_get.side_effect = [fail_response, success_response]
                
                from LigandSimilaritySearcher.sources.lib.similarity import search_similar_compounds
                
                # Should retry and succeed
                result = search_similar_compounds("CCO", n_results=1)
                
                # Verify retry happened
                self.assertEqual(mock_get.call_count, 2)
                self.assertEqual(len(result), 1)


class CSVParsingTests(TestCase):
    """Test CSV parsing and rendering."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_csv_with_all_fields(self) -> None:
        """Test parsing CSV with all expected fields."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        
        csv_content = """CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name
702,CCO,1.0,46.07,-0.31,1,1,20.23,0,-0.77,Ethanol"""
        
        job.out_csv.save("test.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        data = response.context["csv_data"]
        
        self.assertEqual(len(data), 1)
        self.assertEqual(data[0]["CID"], "702")
        self.assertEqual(data[0]["SMILES"], "CCO")
        self.assertEqual(data[0]["Similarity"], "1.0")
        self.assertEqual(data[0]["Name"], "Ethanol")
        
    def test_csv_with_missing_fields(self) -> None:
        """Test graceful handling of CSV with missing fields."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        
        # CSV missing LogS and Name
        csv_content = """CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds
702,CCO,1.0,46.07,-0.31,1,1,20.23,0"""
        
        job.out_csv.save("test.csv", ContentFile(csv_content), save=True)
        
        # Should not crash
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
    def test_csv_with_special_characters(self) -> None:
        """Test CSV with special characters in compound names."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        
        csv_content = """CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name
702,CCO,1.0,46.07,-0.31,1,1,20.23,0,-0.77,"Ethanol, 99%"
1031,CCCO,0.95,60.1,0.25,1,1,20.23,1,-1.01,"1-Propanol (n-Propanol)"
"""
        
        job.out_csv.save("test.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertEqual(response.status_code, 200)
        
        # Verify special characters are handled (CSV parser may convert quotes)
        # Check that names appear in some form
        self.assertContains(response, "Ethanol")
        self.assertContains(response, "Propanol")


class DrugLikenessIndicatorsTests(TestCase):
    """Test drug-likeness indicators in similarity search."""
    
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("tester", password="pass123")
        self.client.force_login(self.user)
        
    def test_optimal_logp_indicator(self) -> None:
        """Test optimal LogP (0-5) shows green indicator."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        
        csv_content = """CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name
702,CCO,1.0,46.07,2.5,1,1,20.23,0,-0.77,Good Drug"""
        
        job.out_csv.save("test.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertContains(response, "bg-success")  # Green badge
        # Check for LogP value
        self.assertContains(response, "2.5")
        
    def test_warning_logp_indicator(self) -> None:
        """Test warning LogP (outside 0-5) shows orange indicator."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        
        csv_content = """CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name
12345,CCCCCCCC,0.90,200.0,6.5,0,0,0.0,7,-5.2,Poor Drug"""
        
        job.out_csv.save("test.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertContains(response, "bg-warning")  # Orange badge
        self.assertContains(response, "High")
        
    def test_perfect_similarity_badge(self) -> None:
        """Test perfect similarity (1.0) shows special badge."""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a="CCO",
            user=self.user
        )
        
        csv_content = """CID,SMILES,Similarity,MW,LogP,HBD,HBA,TPSA,RotBonds,LogS,Name
702,CCO,1.0,46.07,-0.31,1,1,20.23,0,-0.77,Perfect Match"""
        
        job.out_csv.save("test.csv", ContentFile(csv_content), save=True)
        
        response = self.client.get(reverse("chemtools:job_detail", args=[job.pk]))
        self.assertContains(response, "Perfect")
        self.assertContains(response, "bg-success")
