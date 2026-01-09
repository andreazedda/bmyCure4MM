#!/usr/bin/env python3
"""Replace the simulation section in patient_detail.html with decision support version."""

template_path = "clinic/templates/clinic/patient_detail.html"
new_section_path = "clinic/templates/clinic/patient_detail_decision_support.html"

# Read the new section
with open(new_section_path, "r", encoding="utf-8") as f:
    new_section = f.read()

# Read the current template
with open(template_path, "r", encoding="utf-8") as f:
    lines = f.readlines()

# Find the start marker: "{% if latest_simulation_attempt %}"
# Find the end marker:  "{% endif %}" before the next section (patient demographics)

start_line = None
end_line = None
depth = 0

for i, line in enumerate(lines):
    # Look for line 102: "                        {% if latest_simulation_attempt %}"
    if "{% if latest_simulation_attempt %}" in line and i > 90 and i < 110:
        start_line = i
        depth = 1
        continue
    
    if start_line is not None:
        # Track nested if/endif
        if "{% if " in line:
            depth += 1
        elif "{% endif %}" in line:
            depth -= 1
            if depth == 0:
                end_line = i + 1  # Include the endif line
                break

if start_line is None or end_line is None:
    print(f"Could not find section boundaries. Start: {start_line}, End: {end_line}")
    exit(1)

print(f"Replacing lines {start_line + 1} to {end_line} (inclusive)")
print(f"Old section: {end_line - start_line} lines")
print(f"New section: {len(new_section.splitlines())} lines")

# Build new content
new_lines = (
    lines[:start_line] +
    [new_section if not new_section.endswith("\n") else new_section] +
    (["\n"] if not new_section.endswith("\n") else []) +
    lines[end_line:]
)

# Write back
with open(template_path, "w", encoding="utf-8") as f:
    f.writelines(new_lines)

print(f"✅ Template updated successfully!")
print(f"   Backup saved at: {template_path}.backup")
