# Generated migration for enhanced Scenario model with clinical validation fields

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('simulator', '0006_help_articles_expanded'),
    ]

    operations = [
        # Cytogenetics fields
        migrations.AddField(
            model_name='scenario',
            name='del17p',
            field=models.BooleanField(
                default=False,
                help_text='TP53 deletion - high risk',
                verbose_name='del(17p)'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='t_4_14',
            field=models.BooleanField(
                default=False,
                help_text='Translocation (4;14) - high risk',
                verbose_name='t(4;14)'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='t_14_16',
            field=models.BooleanField(
                default=False,
                help_text='Translocation (14;16) - very high risk',
                verbose_name='t(14;16)'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='gain_1q21',
            field=models.BooleanField(
                default=False,
                help_text='1q21 gain/amplification - proliferation advantage',
                verbose_name='1q21 gain'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='hyperdiploid',
            field=models.BooleanField(
                default=False,
                help_text='Hyperdiploid karyotype - standard risk',
                verbose_name='Hyperdiploid'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='t_11_14',
            field=models.BooleanField(
                default=False,
                help_text='Translocation (11;14) - standard risk, better prognosis',
                verbose_name='t(11;14)'
            ),
        ),
        
        # Tumor biology parameters
        migrations.AddField(
            model_name='scenario',
            name='tumor_cell_count',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='Tumor cell count (cells). Range: 1e6-1e12. Typical newly diagnosed: 1e10'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='tumor_growth_rate',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='Tumor growth rate (per day). Range: 0.001-0.1. Typical: 0.01'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='carrying_capacity',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='Maximum tumor burden (cells). Range: 1e11-1e13. Default: 1e12'
            ),
        ),
        
        # Patient characteristics
        migrations.AddField(
            model_name='scenario',
            name='patient_age',
            field=models.IntegerField(
                null=True,
                blank=True,
                help_text='Patient age in years. Range: 18-120'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='ecog_performance_status',
            field=models.IntegerField(
                null=True,
                blank=True,
                choices=[
                    (0, '0 - Fully active'),
                    (1, '1 - Restricted in strenuous activity'),
                    (2, '2 - Ambulatory, self-care'),
                    (3, '3 - Limited self-care'),
                    (4, '4 - Completely disabled'),
                ],
                help_text='ECOG Performance Status (0-4)'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='charlson_comorbidity_index',
            field=models.IntegerField(
                null=True,
                blank=True,
                help_text='Charlson Comorbidity Index. Range: 0-10. Low: 0-1, Moderate: 2-3, High: ≥4'
            ),
        ),
        
        # Laboratory values
        migrations.AddField(
            model_name='scenario',
            name='creatinine_clearance',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='Creatinine clearance (mL/min). Normal: >60, Moderate renal impairment: 30-60, Severe: <30'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='albumin',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='Serum albumin (g/dL). Normal: 3.5-5.0'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='beta2_microglobulin',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='β2-microglobulin (mg/L). Normal: <2, Elevated: 2-5.5, Very high: >5.5',
                verbose_name='β2-microglobulin'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='ldh',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='Lactate dehydrogenase (U/L). Normal: 140-280, Elevated: >280',
                verbose_name='LDH'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='hemoglobin',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='Hemoglobin (g/dL). Normal: M 13-17, F 12-15. Anemia: <10'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='calcium',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='Serum calcium (mg/dL). Normal: 8.5-10.2, Hypercalcemia: >10.5, Severe: >14'
            ),
        ),
        
        # R-ISS staging
        migrations.AddField(
            model_name='scenario',
            name='riss_stage',
            field=models.CharField(
                max_length=16,
                blank=True,
                default='',
                choices=[
                    ('', 'Not specified'),
                    ('I', 'R-ISS I (Low risk)'),
                    ('II', 'R-ISS II (Intermediate risk)'),
                    ('III', 'R-ISS III (High risk)'),
                ],
                help_text='Revised International Staging System stage',
                verbose_name='R-ISS Stage'
            ),
        ),
        
        # Patient archetype from virtual patients
        migrations.AddField(
            model_name='scenario',
            name='patient_archetype',
            field=models.CharField(
                max_length=32,
                blank=True,
                default='',
                choices=[
                    ('', 'Custom/Mixed'),
                    ('nd_standard', 'Newly Diagnosed - Standard Risk'),
                    ('nd_high_risk', 'Newly Diagnosed - High Risk'),
                    ('frail_elderly', 'Frail Elderly (Age >75)'),
                    ('renal_impaired', 'Renal Impairment (CrCl <60)'),
                    ('rr_early', 'Relapsed/Refractory - Early (1-2 lines)'),
                    ('rr_late', 'Relapsed/Refractory - Late (≥3 lines)'),
                    ('aggressive', 'Aggressive/Plasma Cell Leukemia'),
                ],
                help_text='Patient archetype from virtual patient generator'
            ),
        ),
        
        # Calculated difficulty score
        migrations.AddField(
            model_name='scenario',
            name='difficulty_score',
            field=models.FloatField(
                null=True,
                blank=True,
                help_text='Calculated difficulty score (0-100). Auto-populated from difficulty scoring system.'
            ),
        ),
        migrations.AddField(
            model_name='scenario',
            name='difficulty_level',
            field=models.CharField(
                max_length=16,
                blank=True,
                default='',
                choices=[
                    ('', 'Not calculated'),
                    ('easy', 'Easy (0-30)'),
                    ('moderate', 'Moderate (30-50)'),
                    ('hard', 'Hard (50-70)'),
                    ('very_hard', 'Very Hard (70-85)'),
                    ('expert', 'Expert (85-100)'),
                ],
                help_text='Difficulty level category'
            ),
        ),
    ]
