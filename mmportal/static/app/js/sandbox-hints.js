/**
 * 💡 Sandbox Hints System - Contextual Learning Support
 * 
 * Provides smart contextual hints based on user actions to guide learning
 * without being intrusive. Hints appear after specific triggers and can be
 * dismissed but reappear after continued interaction.
 * 
 * @module SandboxHints
 */

(function () {
    'use strict';

    const HINTS_KEY = 'sandbox_hints_dismissed';
    const HINT_TRIGGERS = {
        low_dose: {
            check: () => {
                const lenInput = document.getElementById('id_lenalidomide_dose');
                const borInput = document.getElementById('id_bortezomib_dose');
                if (!lenInput || !borInput) return false;
                
                const lenVal = parseFloat(lenInput.value);
                const borVal = parseFloat(borInput.value);
                
                // Trigger if both doses are very low (less than 30% of typical)
                return lenVal > 0 && lenVal < 10 && borVal > 0 && borVal < 0.5;
            },
            hint: {
                en: "💊 Low doses may not control disease effectively. Try standard doses: Lenalidomide 25mg, Bortezomib 1.3mg/m².",
                it: "💊 Dosi basse potrebbero non controllare la malattia efficacemente. Prova dosi standard: Lenalidomide 25mg, Bortezomib 1.3mg/m²."
            },
            icon: "💊"
        },
        high_healthy_loss: {
            check: () => {
                // Check if previous result had high healthy loss
                const healthyLossCell = document.querySelector('[data-kpi="healthy_loss"]');
                if (!healthyLossCell) return false;
                
                const lossText = healthyLossCell.textContent.trim();
                const lossValue = parseFloat(lossText);
                
                return lossValue > 0.25; // >25% healthy cell loss
            },
            hint: {
                en: "⚠️ High healthy cell loss detected! Consider: reducing doses, adding growth factor support, or increasing rest periods.",
                it: "⚠️ Elevata perdita di cellule sane rilevata! Considera: ridurre le dosi, aggiungere fattori di crescita, o aumentare i periodi di riposo."
            },
            icon: "⚠️"
        },
        no_twin: {
            check: () => {
                const twinCheckbox = document.getElementById('id_enable_digital_twin');
                const attemptCount = document.querySelectorAll('[data-attempt-id]').length;
                
                // Trigger after 3+ simulations without enabling twin
                return twinCheckbox && !twinCheckbox.checked && attemptCount >= 3;
            },
            hint: {
                en: "🎯 Try enabling Digital Twin! It personalizes predictions based on patient biomarkers for more accurate results.",
                it: "🎯 Prova ad abilitare Digital Twin! Personalizza le previsioni basandosi sui biomarker del paziente per risultati più accurati."
            },
            icon: "🎯"
        },
        first_simulation: {
            check: () => {
                const attemptCount = document.querySelectorAll('[data-attempt-id]').length;
                const hasInteractedWithForm = sessionStorage.getItem('form_interacted') === 'true';
                
                // Trigger on first visit or after resetting form
                return attemptCount === 0 && !hasInteractedWithForm;
            },
            hint: {
                en: "🚀 New to bmyCure4MM? Start with a preset regimen like VRd to see the simulator in action!",
                it: "🚀 Nuovo su bmyCure4MM? Inizia con un regime preimpostato come VRd per vedere il simulatore in azione!"
            },
            icon: "🚀"
        },
        extreme_parameters: {
            check: () => {
                const tumorGrowth = document.getElementById('id_tumor_growth_rate');
                const healthyGrowth = document.getElementById('id_healthy_growth_rate');
                
                if (!tumorGrowth || !healthyGrowth) return false;
                
                const tgVal = parseFloat(tumorGrowth.value);
                const hgVal = parseFloat(healthyGrowth.value);
                
                // Trigger if parameters are outside typical biological range
                return (tgVal > 0.05 || hgVal > 0.03);
            },
            hint: {
                en: "🔬 Extreme growth rates detected! Most MM cases have tumor growth 0.015-0.03 and healthy growth 0.01-0.02 day⁻¹.",
                it: "🔬 Tassi di crescita estremi rilevati! La maggior parte dei casi MM ha crescita tumorale 0.015-0.03 e crescita sana 0.01-0.02 day⁻¹."
            },
            icon: "🔬"
        },
        long_horizon: {
            check: () => {
                const horizonInput = document.getElementById('id_time_horizon');
                if (!horizonInput) return false;
                
                const horizon = parseFloat(horizonInput.value);
                return horizon > 240; // >8 months
            },
            hint: {
                en: "⏱️ Long simulation horizons (>8 months) may be less accurate. Consider 180 days for one treatment cycle.",
                it: "⏱️ Orizzonti di simulazione lunghi (>8 mesi) potrebbero essere meno accurati. Considera 180 giorni per un ciclo di trattamento."
            },
            icon: "⏱️"
        }
    };

    class SandboxHints {
        constructor() {
            this.dismissed = this.loadDismissed();
            this.actionCount = 0;
            this.activeHints = [];
            this.panel = null;
            this.button = null;
            
            this.init();
        }

        init() {
            this.createUI();
            this.attachListeners();
            this.checkTriggers();
        }

        createUI() {
            // Create floating hint button
            this.button = document.createElement('button');
            this.button.className = 'sandbox-hint-button';
            this.button.innerHTML = '💡';
            this.button.setAttribute('aria-label', 'Show contextual hints');
            this.button.style.display = 'none';
            this.button.addEventListener('click', () => this.togglePanel());

            // Create hints panel
            this.panel = document.createElement('div');
            this.panel.className = 'sandbox-hints-panel';
            this.panel.innerHTML = `
                <div class="d-flex justify-content-between align-items-center mb-2">
                    <h6 class="mb-0">💡 <span class="t-en">Quick Tips</span><span class="t-it">Suggerimenti Rapidi</span></h6>
                    <button class="btn-close" aria-label="Close"></button>
                </div>
                <div id="hint-items-container"></div>
            `;

            document.body.appendChild(this.button);
            document.body.appendChild(this.panel);

            this.panel.querySelector('.btn-close').addEventListener('click', () => {
                this.hidePanel();
            });
        }

        attachListeners() {
            // Track form interactions
            const form = document.getElementById('simulation-parameters-form');
            if (form) {
                form.addEventListener('input', () => {
                    sessionStorage.setItem('form_interacted', 'true');
                    this.actionCount++;
                    
                    // Check triggers after every 5 actions
                    if (this.actionCount % 5 === 0) {
                        this.checkTriggers();
                    }
                });
            }

            // Check triggers after simulation runs
            document.body.addEventListener('htmx:afterSwap', (event) => {
                if (event.detail.target.id === 'simulation-results') {
                    setTimeout(() => this.checkTriggers(), 1000);
                }
            });

            // Recheck triggers on preset change
            const presetSelect = document.getElementById('id_preset');
            if (presetSelect) {
                presetSelect.addEventListener('change', () => {
                    setTimeout(() => this.checkTriggers(), 500);
                });
            }
        }

        checkTriggers() {
            const newHints = [];

            for (const [key, trigger] of Object.entries(HINT_TRIGGERS)) {
                // Skip if already dismissed
                if (this.dismissed.includes(key)) continue;

                // Check if trigger condition is met
                if (trigger.check()) {
                    const lang = this.getCurrentLang();
                    newHints.push({
                        key: key,
                        icon: trigger.icon,
                        text: trigger.hint[lang] || trigger.hint.en
                    });
                }
            }

            if (newHints.length > 0) {
                this.activeHints = newHints.slice(0, 3); // Max 3 hints at once
                this.renderHints();
                this.showButton();
            } else {
                this.hideButton();
            }
        }

        renderHints() {
            const container = this.panel.querySelector('#hint-items-container');
            container.innerHTML = '';

            this.activeHints.forEach(hint => {
                const hintEl = document.createElement('div');
                hintEl.className = 'hint-item';
                hintEl.innerHTML = `
                    <span class="hint-icon">${hint.icon}</span>
                    <div class="flex-grow-1">
                        <p class="mb-0">${hint.text}</p>
                    </div>
                    <button class="btn btn-sm btn-link text-muted p-0" 
                            data-hint-key="${hint.key}"
                            aria-label="Dismiss hint">×</button>
                `;

                const dismissBtn = hintEl.querySelector('button');
                dismissBtn.addEventListener('click', () => {
                    this.dismissHint(hint.key);
                    hintEl.remove();
                    
                    // Hide panel if no hints left
                    if (this.panel.querySelectorAll('.hint-item').length === 0) {
                        this.hidePanel();
                        this.hideButton();
                    }
                });

                container.appendChild(hintEl);
            });
        }

        dismissHint(key) {
            if (!this.dismissed.includes(key)) {
                this.dismissed.push(key);
                this.saveDismissed();
            }
        }

        showButton() {
            this.button.style.display = 'flex';
        }

        hideButton() {
            this.button.style.display = 'none';
        }

        togglePanel() {
            if (this.panel.classList.contains('show')) {
                this.hidePanel();
            } else {
                this.showPanel();
            }
        }

        showPanel() {
            this.panel.classList.add('show');
        }

        hidePanel() {
            this.panel.classList.remove('show');
        }

        loadDismissed() {
            try {
                const stored = localStorage.getItem(HINTS_KEY);
                return stored ? JSON.parse(stored) : [];
            } catch (e) {
                return [];
            }
        }

        saveDismissed() {
            try {
                localStorage.setItem(HINTS_KEY, JSON.stringify(this.dismissed));
            } catch (e) {
                console.warn('Failed to save dismissed hints', e);
            }
        }

        getCurrentLang() {
            return (window.getCurrentLang && window.getCurrentLang()) || 'en';
        }

        // Public API: Reset all dismissed hints (for testing/admin)
        resetDismissed() {
            this.dismissed = [];
            this.saveDismissed();
            this.checkTriggers();
        }
    }

    // Initialize when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', () => {
            window.SandboxHints = new SandboxHints();
        });
    } else {
        window.SandboxHints = new SandboxHints();
    }

    // Expose reset function globally for debugging
    window.resetSandboxHints = function() {
        if (window.SandboxHints) {
            window.SandboxHints.resetDismissed();
            console.log('Sandbox hints reset');
        }
    };
})();
