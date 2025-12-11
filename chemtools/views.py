from __future__ import annotations

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.http import HttpRequest, HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.views.decorators.http import require_GET

from . import forms, models
from .tasks import (
    run_binding_viz_job,
    run_drug_params_job,
    run_similarity_job,
)
from .pdb_api_client import enrich_pdb_metadata_for_view


def _enqueue(task, job: models.ChemJob, *params) -> tuple[bool, bool]:
    """Attempt to queue a Celery task; fallback to synchronous execution.

    Returns a tuple ``(queued, failed)`` where ``queued`` indicates the task
    was handed to Celery, and ``failed`` reflects the synchronous run status.
    """
    try:
        task.delay(job.pk, *params)
        return True, False
    except Exception:  # pragma: no cover - fallback when broker unavailable
        task.run(job.pk, *params)
        job.refresh_from_db()
        failed = bool(job.log and job.log.startswith("ERROR"))
        return False, failed


@login_required
def tools_home(request: HttpRequest) -> HttpResponse:
    # Filter jobs by current user for privacy
    jobs = models.ChemJob.objects.filter(user=request.user).select_related("user").order_by('-created')[:50]
    return render(
        request,
        "chemtools/tools_home.html",
        {
            "jobs": jobs,
        },
    )


@login_required
def drug_params(request: HttpRequest) -> HttpResponse:
    form = forms.DrugParamForm(request.POST or None)
    if request.method == "POST" and form.is_valid():
        smiles = (form.cleaned_data.get("smiles") or "").strip()
        cid_value = form.cleaned_data.get("cid")
        cid_str = str(cid_value) if cid_value else ""
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.PARAM,
            input_a=smiles,
            input_b=cid_str,
            user=request.user if request.user.is_authenticated else None,
        )
        queued, failed = _enqueue(run_drug_params_job, job, smiles, cid_str)
        if queued:
            messages.info(request, "Job queued. Refresh this page to monitor progress.")
        elif failed:
            messages.error(request, "Drug parameter evaluation failed. See log for details.")
        else:
            messages.success(request, "Drug parameter evaluation completed.")
        return redirect(reverse("chemtools:tools_home"))
    return render(request, "chemtools/drug_params.html", {"form": form})


@login_required
def binding_viz(request: HttpRequest) -> HttpResponse:
    import logging
    logger = logging.getLogger(__name__)
    
    logger.info(f"[DEBUG] 🌐 binding_viz view called - Method: {request.method}, User: {request.user}")
    
    form = forms.BindingVizForm(request.POST or None)
    if request.method == "POST" and form.is_valid():
        pdb_id = form.cleaned_data["pdb_id"]
        ligand = form.cleaned_data.get("ligand", "")
        
        logger.info(f"[DEBUG] 📝 Form validated - PDB ID: {pdb_id}, Ligand: {ligand or '(none)'}")
        
        # Capture API preferences
        # BooleanField with required=False returns False when unchecked
        # Use the cleaned_data directly as it reflects user's actual choices
        api_prefs = {
            'fetch_validation': form.cleaned_data.get('fetch_validation', True),
            'fetch_interactions': form.cleaned_data.get('fetch_interactions', True),
            'fetch_drug_info': form.cleaned_data.get('fetch_drug_info', True),
            'fetch_clinical_trials': form.cleaned_data.get('fetch_clinical_trials', True),
            'fetch_publications': form.cleaned_data.get('fetch_publications', True),
            'fetch_protein_network': form.cleaned_data.get('fetch_protein_network', True),
            'fetch_pathways': form.cleaned_data.get('fetch_pathways', True),
            'mm_focus_resistance': form.cleaned_data.get('mm_focus_resistance', False),
            'mm_focus_combinations': form.cleaned_data.get('mm_focus_combinations', False),
            'mm_focus_toxicity': form.cleaned_data.get('mm_focus_toxicity', False),
        }
        
        logger.debug(f"[DEBUG] 📊 API Preferences captured: {api_prefs}")
        enabled_apis = [k for k, v in api_prefs.items() if v]
        logger.info(f"[DEBUG] ✅ Enabled APIs ({len(enabled_apis)}): {', '.join(enabled_apis)}")
        
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.BIND,
            input_a=pdb_id,
            input_b=ligand,
            api_preferences=api_prefs,
            user=request.user if request.user.is_authenticated else None,
        )
        
        logger.info(f"[DEBUG] 💾 Job created - ID: {job.pk}, Kind: BIND")
        logger.debug(f"[DEBUG] 🔗 Job redirect URL: {reverse('chemtools:job_detail', args=[job.pk])}")
        
        logger.info(f"[DEBUG] 🚀 Attempting to enqueue job {job.pk}...")
        queued, failed = _enqueue(run_binding_viz_job, job, pdb_id, ligand)
        
        if queued:
            logger.info(f"[DEBUG] ✅ Job {job.pk} queued successfully (Celery)")
            messages.info(request, "Binding visualization queued. Refresh shortly for results.")
        elif failed:
            logger.error(f"[DEBUG] ❌ Job {job.pk} failed during synchronous execution")
            messages.error(request, "Binding visualization failed. Check the job log for details.")
        else:
            logger.info(f"[DEBUG] ✅ Job {job.pk} completed synchronously")
            messages.success(request, "Binding visualization completed.")
        
        return redirect(reverse("chemtools:tools_home"))
    
    logger.debug(f"[DEBUG] 📄 Rendering binding_viz form (GET or invalid POST)")
    return render(request, "chemtools/binding_viz.html", {"form": form})


@login_required
def similarity(request: HttpRequest) -> HttpResponse:
    form = forms.SimilarityForm(request.POST or None)
    if request.method == "POST" and form.is_valid():
        smiles = form.cleaned_data["smiles"].strip()
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a=smiles,
            user=request.user if request.user.is_authenticated else None,
        )
        queued, failed = _enqueue(run_similarity_job, job, smiles)
        if queued:
            messages.info(request, "Similarity search queued. Refresh shortly for results.")
        elif failed:
            messages.error(request, "Similarity search failed. Review the job log for details.")
        else:
            messages.success(request, "Similarity search completed.")
        return redirect(reverse("chemtools:tools_home"))
    return render(request, "chemtools/similarity.html", {"form": form})


@login_required
def retry_job(request: HttpRequest, pk: int) -> HttpResponse:
    job = get_object_or_404(models.ChemJob, pk=pk)
    if job.kind == models.ChemJob.PARAM:
        queued, failed = _enqueue(
            run_drug_params_job,
            job,
            job.input_a or "",
            job.input_b or "",
        )
    elif job.kind == models.ChemJob.BIND:
        queued, failed = _enqueue(
            run_binding_viz_job,
            job,
            job.input_a,
            job.input_b or "",
        )
    else:
        queued, failed = _enqueue(run_similarity_job, job, job.input_a)

    if queued:
        messages.info(request, "Job retry queued. Refresh to monitor progress.")
    elif failed:
        messages.error(request, "Job retry failed immediately. Check the log for details.")
    else:
        messages.success(request, "Job retry completed successfully.")
    return redirect(reverse("chemtools:tools_home"))


@login_required
@require_GET
def job_status(request: HttpRequest, pk: int) -> JsonResponse:
    job = get_object_or_404(models.ChemJob, pk=pk)
    status_text, status_variant = job.status_label()
    data = {
        "status": status_text,
        "variant": status_variant,
        "has_html": bool(job.out_html),
        "has_csv": bool(job.out_csv),
        "html_url": job.out_html.url if job.out_html else None,
        "csv_url": job.out_csv.url if job.out_csv else None,
        "thumbnail_url": job.thumbnail_url(),
        "progress_percent": job.progress_percent,
        "progress_message": job.progress_message,
    }
    return JsonResponse(data)


@login_required
def job_detail(request: HttpRequest, pk: int) -> HttpResponse:
    """Integrated results viewer for all job types with enhanced analytics."""
    import logging
    import time
    logger = logging.getLogger(__name__)
    
    start_time = time.time()
    logger.info(f"[DEBUG] 📋 job_detail view called - Job ID: {pk}, User: {request.user}")
    
    job = get_object_or_404(models.ChemJob, pk=pk, user=request.user)
    logger.info(f"[DEBUG] 💼 Job loaded - Kind: {job.kind}, Created: {job.created}")
    logger.debug(f"[DEBUG]   Input A: {job.input_a}, Input B: {job.input_b}")
    logger.debug(f"[DEBUG]   Has HTML: {bool(job.out_html)}, Has CSV: {bool(job.out_csv)}")
    logger.debug(f"[DEBUG]   API Preferences: {job.api_preferences if hasattr(job, 'api_preferences') else 'N/A'}")
    
    # Load content based on job type
    html_content = None
    csv_data = None
    pdb_metadata = None
    binding_analysis = None
    
    if job.out_html:
        logger.debug(f"[DEBUG] 📄 Loading HTML output...")
        try:
            with job.out_html.open('r') as f:
                html_content = f.read()
            logger.info(f"[DEBUG] ✅ HTML loaded - Size: {len(html_content)} chars, Lines: {html_content.count(chr(10))}")
        except Exception as e:
            logger.error(f"[DEBUG] ❌ Failed to load HTML: {e}")
            html_content = None
    
    if job.out_csv:
        logger.debug(f"[DEBUG] 📊 Loading CSV output...")
        try:
            import csv as csv_module
            with job.out_csv.open('r') as f:
                reader = csv_module.DictReader(f)
                csv_data = list(reader)
            logger.info(f"[DEBUG] ✅ CSV loaded - Rows: {len(csv_data)}, Columns: {list(csv_data[0].keys()) if csv_data else []}")
        except Exception as e:
            logger.error(f"[DEBUG] ❌ Failed to load CSV: {e}")
            csv_data = None
    
    # Enhanced metadata for binding jobs - NOW WITH API INTEGRATION
    if job.kind == models.ChemJob.BIND and job.input_a:
        pdb_id = job.input_a.upper()
        logger.info(f"[DEBUG] 🔬 Processing binding job metadata for PDB: {pdb_id}")
        
        # Basic PDB metadata URLs
        pdb_metadata = {
            'pdb_id': pdb_id,
            'rcsb_url': f'https://www.rcsb.org/structure/{pdb_id}',
            'pdbe_url': f'https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_id}',
            'pdb_redo_url': f'https://pdb-redo.eu/db/{pdb_id}',
            'alphafold_similar': f'https://alphafold.ebi.ac.uk/search/text/{pdb_id}',
        }
        logger.debug(f"[DEBUG]   Basic URLs configured: {list(pdb_metadata.keys())}")
        
        # Get API preferences from job (use defaults if not set)
        api_prefs = job.api_preferences if hasattr(job, 'api_preferences') and job.api_preferences else None
        if api_prefs:
            enabled_count = sum(1 for v in api_prefs.values() if v)
            logger.info(f"[DEBUG] 🔧 API preferences loaded: {enabled_count} sources enabled")
            logger.debug(f"[DEBUG]   Preferences: {api_prefs}")
        else:
            logger.warning(f"[DEBUG] ⚠️ No API preferences found, using defaults")
        
        # Enrich with API data (includes PDB summary, ligands, validation, etc.)
        logger.info(f"[DEBUG] 🌐 Enriching metadata with API data...")
        enrich_start = time.time()
        api_enriched_data = enrich_pdb_metadata_for_view(pdb_id, ligand_id=job.input_b, api_prefs=api_prefs)
        enrich_time = time.time() - enrich_start
        logger.info(f"[DEBUG] ✅ Metadata enrichment completed in {enrich_time:.2f}s")
        if api_enriched_data:
            logger.debug(f"[DEBUG]   Enriched data keys: {list(api_enriched_data.keys())}")
        else:
            logger.warning(f"[DEBUG] ⚠️ No enriched data returned")
        
        # Merge API data into pdb_metadata
        if api_enriched_data:
            pdb_metadata.update(api_enriched_data)
            
            # Extract binding_analysis from API data for template compatibility
            if 'binding_analysis' in api_enriched_data:
                binding_analysis = api_enriched_data['binding_analysis']
                logger.debug(f"[DEBUG]   Binding analysis extracted: {list(binding_analysis.keys()) if binding_analysis else 'empty'}")
    
    elapsed = time.time() - start_time
    logger.info(f"[DEBUG] ⏱️ job_detail view completed in {elapsed:.2f}s")
    logger.info(f"[DEBUG] 📦 Context prepared with: html_content={bool(html_content)}, csv_data={bool(csv_data)}, pdb_metadata={bool(pdb_metadata)}, binding_analysis={bool(binding_analysis)}")
    
    return render(
        request,
        "chemtools/job_detail.html",
        {
            "job": job,
            "html_content": html_content,
            "csv_data": csv_data,
            "pdb_metadata": pdb_metadata,
            "binding_analysis": binding_analysis,
        },
    )
