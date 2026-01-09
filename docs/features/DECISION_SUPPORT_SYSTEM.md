# 🎯 AI-Powered Decision Support System

> Intelligent treatment optimization assistant for multiple myeloma therapy planning

## Overview

The Decision Support System is an AI-powered assistant that automatically analyzes PKPD simulation results and provides actionable, priority-ranked recommendations to optimize treatment outcomes. It transforms raw simulation data into clear, implementable clinical insights.

## Key Features

### 🤖 Automatic Analysis
- Real-time interpretation of simulation results
- Smart detection of efficacy vs toxicity balance
- Pattern recognition across multiple metrics (tumor reduction, healthy cell loss, time to recurrence)

### 📊 Priority-Based Recommendations
Four priority levels guide attention to critical issues:
- 🚨 **Critical**: Immediate action required (e.g., tumor growth, severe toxicity)
- ⚠️ **High**: Important optimization needed (e.g., high toxicity, poor efficacy)
- ⚙️ **Medium**: Suggested improvements (e.g., early recurrence risk)
- ✅ **Low**: Fine-tuning opportunities (e.g., already favorable results)

### 💡 Actionable Guidance
Each recommendation includes:
- **Issue Description**: What problem was detected
- **Specific Action**: Numerical guidance (e.g., "Reduce doses by 20-30%", "Extend horizon to 224 days")
- **Rationale**: Why this action is recommended
- **Direct Link**: One-click access to the scenario for implementation

### 🎨 Smart UI Behavior
- **Auto-Opening Accordion**: Recommendations panel opens automatically when warnings detected
- **Color-Coded Alerts**: Green (favorable) / Yellow (caution) / Red (unfavorable)
- **Visual Priority**: Critical items highlighted with danger background
- **Badge System**: Quick visual identification of metric quality

## How It Works

### 1. Simulation Execution
When you run a PKPD simulation in the Simulator module, the system:
1. Executes the mathematical model with your chosen parameters
2. Calculates key performance indicators:
   - Tumor reduction (efficacy)
   - Healthy cell loss (toxicity)
   - Time to recurrence (durability)

### 2. Intelligent Analysis
The `_interpret_latest_simulation()` function analyzes results using:

**Efficacy Thresholds:**
- ≥ 0.50: Good tumor reduction
- 0.30-0.50: Moderate response
- < 0.30: Poor efficacy
- < 0: Tumor growth (critical)

**Toxicity Thresholds:**
- < 0.20: Acceptable healthy cell loss
- 0.20-0.30: Borderline toxicity
- ≥ 0.30: High toxicity (concerning)

**Durability Thresholds:**
- ≥ 180 days: Good durability
- 90-180 days: Moderate
- < 90 days: Early recurrence risk

### 3. Recommendation Generation

The system generates specific recommendations based on detected patterns:

#### Scenario 1: High Toxicity (≥30% healthy loss)
```
🚨 Priority: HIGH/CRITICAL
Issue: "High toxicity (healthy cell loss ≥30%)"
Action: "Reduce drug doses by 20-30% or shorten time horizon"
Rationale: "Too much damage to healthy plasma cells. Lower doses preserve immune function."
```

#### Scenario 2: Moderate Toxicity (20-30% healthy loss)
```
⚠️ Priority: MEDIUM
Issue: "Moderate toxicity (healthy cell loss 20-30%)"
Action: "Consider reducing doses by 10-15% if patient shows clinical toxicity signs"
Rationale: "Borderline toxicity. Monitor closely; reduce if side effects appear."
```

#### Scenario 3: Poor Efficacy (<30% tumor reduction)
```
📈 Priority: HIGH
Issue: "Low tumor reduction (<30%)"
Action: "Increase drug doses by 15-25% or extend time horizon to 224-280 days"
Rationale: "More aggressive therapy may improve response if toxicity remains acceptable."
```

#### Scenario 4: Tumor Growth (negative reduction)
```
🚨 Priority: CRITICAL
Issue: "Tumor growth (negative reduction)"
Action: "Switch to a different regimen or significantly increase doses"
Rationale: "Current regimen is ineffective. Consider alternative drug combinations."
```

#### Scenario 5: Early Recurrence (<90 days)
```
⏱️ Priority: MEDIUM
Issue: "Early recurrence predicted (<90 days)"
Action: "Extend time horizon to 224-280 days to simulate longer treatment"
Rationale: "Longer therapy duration may delay recurrence and improve durability."
```

#### Scenario 6: Favorable Balance
```
✅ Priority: LOW
Issue: "Favorable balance (good efficacy, acceptable toxicity)"
Action: "Fine-tune by testing ±10% dose variations or compare alternative regimens"
Rationale: "Current settings look promising. Minor adjustments may further optimize."
```

### 4. User Interface Integration

#### Quick Read Banner
At the top of the patient page, a color-coded banner provides instant assessment:
- 🎯 **Green**: "Favorable overall — efficacy and toxicity look balanced"
- ⚠️ **Yellow**: "Mixed signal — review recommendations below"
- 🚨 **Red**: "Unfavorable signal — see action plan below"

#### Recommendations Accordion
- **Automatic Opening**: Opens immediately if caution/bad status detected
- **Priority Sorting**: Critical items appear first
- **Visual Hierarchy**: Critical items have red background, high items have yellow
- **Direct Action Buttons**: Each recommendation has "🔧 Go to Scenario & Implement" button

#### Metrics Table
Each key metric displayed with:
- Numerical value (e.g., 0.645 tumor reduction)
- Quality badge (good/caution/bad)
- Interpretation guide ("↑ Higher is better", "↓ Lower is better")
- Threshold explanations

## Implementation Guide

### For Users

**Step 1: Run Simulation**
1. Go to patient page
2. Click "Start simulation now" → Opens Simulator
3. Select/create scenario with desired regimen and doses
4. Click "Run simulation" in Simulation Panel
5. Wait 10-30 seconds for results

**Step 2: Review Recommendations**
1. Return to patient page
2. Check Quick Read banner (green/yellow/red)
3. If warnings present, accordion opens automatically showing recommendations
4. Read priority-ordered action items

**Step 3: Implement Changes**
1. Click "🔧 Go to Scenario & Implement" button on any recommendation
2. Opens the scenario page directly
3. Modify drug doses or time horizon as suggested
4. Re-run simulation
5. Return to patient page to see updated results

**Step 4: Iterate**
- Compare before/after results
- Continue optimizing until favorable balance achieved
- Save optimal scenario for reference

### For Developers

**Backend Function: `_interpret_latest_simulation()`**

Location: `clinic/views.py`

```python
def _interpret_latest_simulation(
    summary: dict[str, object] | None, 
    parameters: dict[str, object] | None
) -> dict[str, object] | None:
    """
    Comprehensive decision support: interpretation + actionable recommendations.
    
    Args:
        summary: Simulation results with tumor_reduction, healthy_loss, time_to_recurrence
        parameters: Simulation parameters including time_horizon
    
    Returns:
        {
            'overall': 'good'|'caution'|'bad',
            'tumor_reduction_label': 'good'|'caution'|'bad',
            'healthy_loss_label': 'good'|'caution'|'bad',
            'time_to_recurrence_label': 'good'|'caution'|'bad',
            'recommendations': [
                {
                    'priority': 'critical'|'high'|'medium'|'low',
                    'issue_en': 'English description',
                    'issue_it': 'Italian description',
                    'action_en': 'English action',
                    'action_it': 'Italian action',
                    'rationale_en': 'English rationale',
                    'rationale_it': 'Italian rationale',
                    'icon': '🚨'
                },
                ...
            ],
            'has_recommendations': True|False
        }
    """
```

**Template Integration: `patient_detail.html`**

Key template variables:
- `latest_simulation_interpretation`: Output of `_interpret_latest_simulation()`
- `latest_simulation_scenario_url`: Direct link to scenario for implementation
- `latest_simulation_summary`: Raw metrics from simulation

**Rendering Logic:**
1. Check `latest_simulation_interpretation.overall`
2. Render color-coded banner
3. If `has_recommendations=True`, render accordion
4. Set `aria-expanded="true"` and class `show` if overall != 'good'
5. For each recommendation, render card with priority styling
6. Include direct link button with `latest_simulation_scenario_url`

## Best Practices

### Clinical Use
1. **Always consider patient context**: Recommendations are heuristic guides, not clinical prescriptions
2. **Monitor real patient response**: Model predictions are approximations
3. **Balance efficacy vs toxicity**: Prioritize patient safety and quality of life
4. **Use for scenario exploration**: Compare multiple strategies before decision
5. **Document rationale**: Record why specific doses/regimens were chosen

### Technical Configuration

**Threshold Customization:**
If you need to adjust thresholds for your institution, modify in `clinic/views.py`:

```python
# Efficacy thresholds
tumor_label = band_label(tumor_reduction, good_ge=0.50, bad_lt=0.0)

# Toxicity thresholds
if healthy_loss < 0.20:
    healthy_label = "good"
elif healthy_loss < 0.30:
    healthy_label = "caution"
else:
    healthy_label = "bad"

# Durability thresholds
if time_to_recurrence >= 180:
    recurrence_label = "good"
elif time_to_recurrence >= 90:
    recurrence_label = "caution"
else:
    recurrence_label = "bad"
```

**Adding New Scenarios:**
To add custom recommendation scenarios:

```python
# In _interpret_latest_simulation() function
if <your_condition>:
    recommendations.append({
        "issue_en": "English issue description",
        "issue_it": "Italian issue description",
        "action_en": "English action with specific values",
        "action_it": "Italian action with specific values",
        "rationale_en": "English rationale",
        "rationale_it": "Italian rationale",
        "icon": "🔬",  # Choose appropriate emoji
        "priority": "high",  # critical/high/medium/low
    })
```

## Troubleshooting

### Recommendations Not Appearing
1. Verify simulation completed successfully
2. Check `latest_simulation_summary` contains expected keys
3. Ensure `latest_simulation_interpretation` is passed to template
4. Verify `has_recommendations=True` in interpretation dict

### Accordion Not Auto-Opening on Warning
1. Check `overall` value is 'caution' or 'bad'
2. Verify template has correct conditional: `{% if latest_simulation_interpretation.overall != 'good' %}show{% endif %}`
3. Ensure Bootstrap JS is loaded

### "Go to Scenario" Button Not Working
1. Verify `latest_simulation_scenario_url` is set in view context
2. Check scenario still exists in database
3. Ensure user has permission to view scenario

### Wrong Priority or Missing Scenarios
1. Review threshold values in `_interpret_latest_simulation()`
2. Check all expected scenarios are implemented
3. Verify priority sorting logic: `priority_order = {"critical": 0, "high": 1, "medium": 2, "low": 3}`

## API Integration

For programmatic access to recommendations:

```python
from clinic.views import _interpret_latest_simulation

summary = {
    'tumor_reduction': 0.35,
    'healthy_loss': 0.28,
    'time_to_recurrence': 120
}

parameters = {
    'time_horizon': 168
}

interpretation = _interpret_latest_simulation(summary, parameters)

print(interpretation['overall'])  # 'caution'
print(len(interpretation['recommendations']))  # Number of recommendations
for rec in interpretation['recommendations']:
    print(f"{rec['priority']}: {rec['action_en']}")
```

## Future Enhancements

Planned features:
- [ ] Machine learning-based threshold adaptation
- [ ] Historical outcome tracking for recommendation validation
- [ ] Multi-objective optimization suggestions
- [ ] Drug interaction warnings
- [ ] Patient-specific risk factor integration
- [ ] Export recommendations as PDF reports
- [ ] Real-time collaboration notes on recommendations

## References

- [PKPD Models Documentation](MATHEMATICAL_MODELS_DOCUMENTATION.md)
- [Simulator User Guide](../en/simulator.md)
- [Patient Twin System](../en/patient_twin.md)
- [Clinical Decision Support Best Practices](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6372467/)

---

**Last Updated:** January 2026  
**Version:** 2.0  
**Contact:** For questions or feedback, open an issue on GitHub
