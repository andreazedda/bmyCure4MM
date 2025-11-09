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
    form = forms.BindingVizForm(request.POST or None)
    if request.method == "POST" and form.is_valid():
        pdb_id = form.cleaned_data["pdb_id"]
        ligand = form.cleaned_data.get("ligand", "")
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.BIND,
            input_a=pdb_id,
            input_b=ligand,
            user=request.user if request.user.is_authenticated else None,
        )
        queued, failed = _enqueue(run_binding_viz_job, job, pdb_id, ligand)
        if queued:
            messages.info(request, "Binding visualization queued. Refresh shortly for results.")
        elif failed:
            messages.error(request, "Binding visualization failed. Check the job log for details.")
        else:
            messages.success(request, "Binding visualization completed.")
        return redirect(reverse("chemtools:tools_home"))
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
    """Integrated results viewer for all job types."""
    job = get_object_or_404(models.ChemJob, pk=pk, user=request.user)
    
    # Load content based on job type
    html_content = None
    csv_data = None
    
    if job.out_html:
        try:
            with job.out_html.open('r') as f:
                html_content = f.read()
        except Exception:
            html_content = None
    
    if job.out_csv:
        try:
            import csv as csv_module
            with job.out_csv.open('r') as f:
                reader = csv_module.DictReader(f)
                csv_data = list(reader)
        except Exception:
            csv_data = None
    
    return render(
        request,
        "chemtools/job_detail.html",
        {
            "job": job,
            "html_content": html_content,
            "csv_data": csv_data,
        },
    )
