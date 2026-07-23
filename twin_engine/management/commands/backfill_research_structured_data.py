from __future__ import annotations

from twin_engine.management.commands.import_research_dataset import (
    Command as ImportResearchDatasetCommand,
)


class Command(ImportResearchDatasetCommand):
    """
    Compatibility alias for the former hard-coded research-data backfill.

    Structured clinical research records must now be supplied explicitly
    through a versioned private dataset and manifest.
    """

    help = (
        "Compatibility alias for import_research_dataset. "
        "Requires an explicit versioned private dataset; "
        "no clinical records are embedded in application code."
    )
