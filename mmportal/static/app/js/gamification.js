/**
 * Gamification System for MM Portal
 * Tracks user progress, achievements, and provides educational guidance
 */

(function (root) {
  'use strict';

  // 🎮 User Progress Storage (localStorage)
  const STORAGE_KEY = 'mmportal_progress';
  
  // 🏆 Achievement Definitions
  const ACHIEVEMENTS = {
    first_simulation: {
      id: 'first_simulation',
      icon: '🎯',
      title: { en: 'First Steps', it: 'Primi Passi' },
      description: { en: 'Completed your first simulation!', it: 'Completata la prima simulazione!' },
      points: 10,
    },
    explorer: {
      id: 'explorer',
      icon: '🔍',
      title: { en: 'Explorer', it: 'Esploratore' },
      description: { en: 'Tried 5 different drug combinations', it: 'Provate 5 diverse combinazioni di farmaci' },
      points: 25,
    },
    safety_first: {
      id: 'safety_first',
      icon: '🛡️',
      title: { en: 'Safety First', it: 'Sicurezza Prima' },
      description: { en: 'Maintained healthy cells above 70%', it: 'Mantenute cellule sane sopra il 70%' },
      points: 30,
    },
    tumor_buster: {
      id: 'tumor_buster',
      icon: '💥',
      title: { en: 'Tumor Buster', it: 'Distruttore di Tumori' },
      description: { en: 'Achieved >90% tumor reduction', it: 'Raggiunta riduzione tumorale >90%' },
      points: 50,
    },
    perfect_balance: {
      id: 'perfect_balance',
      icon: '⚖️',
      title: { en: 'Perfect Balance', it: 'Equilibrio Perfetto' },
      description: { en: '>80% tumor reduction + >60% healthy cells', it: '>80% riduzione tumorale + >60% cellule sane' },
      points: 75,
    },
    knowledge_seeker: {
      id: 'knowledge_seeker',
      icon: '📚',
      title: { en: 'Knowledge Seeker', it: 'Cercatore di Conoscenza' },
      description: { en: 'Read 10 help articles', it: 'Letti 10 articoli di aiuto' },
      points: 20,
    },
    advanced_researcher: {
      id: 'advanced_researcher',
      icon: '🔬',
      title: { en: 'Advanced Researcher', it: 'Ricercatore Avanzato' },
      description: { en: 'Completed 20 simulations', it: 'Completate 20 simulazioni' },
      points: 100,
    },
  };

  // 🎓 Skill Levels
  const SKILL_LEVELS = [
    { level: 1, name: { en: 'Novice', it: 'Novizio' }, minPoints: 0, icon: '🌱' },
    { level: 2, name: { en: 'Student', it: 'Studente' }, minPoints: 50, icon: '📖' },
    { level: 3, name: { en: 'Researcher', it: 'Ricercatore' }, minPoints: 150, icon: '🔬' },
    { level: 4, name: { en: 'Expert', it: 'Esperto' }, minPoints: 300, icon: '🎓' },
    { level: 5, name: { en: 'Master', it: 'Maestro' }, minPoints: 500, icon: '👨‍⚕️' },
    { level: 6, name: { en: 'Legend', it: 'Leggenda' }, minPoints: 1000, icon: '🏆' },
  ];

  // 📊 Progress Tracker
  class ProgressTracker {
    constructor() {
      this.load();
    }

    load() {
      const stored = localStorage.getItem(STORAGE_KEY);
      if (stored) {
        try {
          this.data = JSON.parse(stored);
        } catch (e) {
          this.data = this.getDefaultData();
        }
      } else {
        this.data = this.getDefaultData();
      }
    }

    getDefaultData() {
      return {
        points: 0,
        achievements: [],
        simulations_completed: 0,
        help_articles_read: [],
        drug_combinations_tried: [],
        best_tumor_reduction: 0,
        best_healthy_preservation: 0,
        first_visit: new Date().toISOString(),
        last_active: new Date().toISOString(),
      };
    }

    save() {
      this.data.last_active = new Date().toISOString();
      localStorage.setItem(STORAGE_KEY, JSON.stringify(this.data));
    }

    // 🎯 Track Simulation Completion
    recordSimulation(results) {
      this.data.simulations_completed += 1;
      
      // Track best results
      if (results.tumor_reduction > this.data.best_tumor_reduction) {
        this.data.best_tumor_reduction = results.tumor_reduction;
      }
      if (results.healthy_cells > this.data.best_healthy_preservation) {
        this.data.best_healthy_preservation = results.healthy_cells;
      }

      // Track drug combinations
      const combo = [results.drugs].sort().join('+');
      if (!this.data.drug_combinations_tried.includes(combo)) {
        this.data.drug_combinations_tried.push(combo);
      }

      this.checkAchievements(results);
      this.save();
    }

    // 📚 Track Help Article Read
    recordHelpRead(slug) {
      if (!this.data.help_articles_read.includes(slug)) {
        this.data.help_articles_read.push(slug);
        this.data.points += 2; // Small reward for learning
        this.checkAchievements();
        this.save();
      }
    }

    // 🏆 Check and Unlock Achievements
    checkAchievements(results) {
      const newAchievements = [];

      // First simulation
      if (this.data.simulations_completed === 1 && !this.hasAchievement('first_simulation')) {
        newAchievements.push(this.unlockAchievement('first_simulation'));
      }

      // Explorer
      if (this.data.drug_combinations_tried.length >= 5 && !this.hasAchievement('explorer')) {
        newAchievements.push(this.unlockAchievement('explorer'));
      }

      // Advanced researcher
      if (this.data.simulations_completed >= 20 && !this.hasAchievement('advanced_researcher')) {
        newAchievements.push(this.unlockAchievement('advanced_researcher'));
      }

      // Knowledge seeker
      if (this.data.help_articles_read.length >= 10 && !this.hasAchievement('knowledge_seeker')) {
        newAchievements.push(this.unlockAchievement('knowledge_seeker'));
      }

      // Result-based achievements
      if (results) {
        if (results.healthy_cells > 0.7 && !this.hasAchievement('safety_first')) {
          newAchievements.push(this.unlockAchievement('safety_first'));
        }
        if (results.tumor_reduction > 0.9 && !this.hasAchievement('tumor_buster')) {
          newAchievements.push(this.unlockAchievement('tumor_buster'));
        }
        if (results.tumor_reduction > 0.8 && results.healthy_cells > 0.6 && !this.hasAchievement('perfect_balance')) {
          newAchievements.push(this.unlockAchievement('perfect_balance'));
        }
      }

      // Show achievement notifications
      newAchievements.forEach(function (achievement) {
        window.showAchievement(achievement);
      });

      return newAchievements;
    }

    unlockAchievement(id) {
      const achievement = ACHIEVEMENTS[id];
      if (achievement && !this.data.achievements.includes(id)) {
        this.data.achievements.push(id);
        this.data.points += achievement.points;
        return achievement;
      }
      return null;
    }

    hasAchievement(id) {
      return this.data.achievements.includes(id);
    }

    getLevel() {
      const points = this.data.points;
      for (var i = SKILL_LEVELS.length - 1; i >= 0; i--) {
        if (points >= SKILL_LEVELS[i].minPoints) {
          return SKILL_LEVELS[i];
        }
      }
      return SKILL_LEVELS[0];
    }

    getNextLevel() {
      const currentLevel = this.getLevel();
      const currentIndex = SKILL_LEVELS.indexOf(currentLevel);
      if (currentIndex < SKILL_LEVELS.length - 1) {
        return SKILL_LEVELS[currentIndex + 1];
      }
      return null;
    }

    getProgress() {
      const current = this.getLevel();
      const next = this.getNextLevel();
      if (!next) {
        return { percentage: 100, pointsToNext: 0 };
      }
      const pointsInLevel = this.data.points - current.minPoints;
      const pointsNeeded = next.minPoints - current.minPoints;
      return {
        percentage: Math.round((pointsInLevel / pointsNeeded) * 100),
        pointsToNext: next.minPoints - this.data.points,
      };
    }

    // 🎁 Get Personalized Tips
    getTip() {
      const lang = window.getCurrentLang ? window.getCurrentLang() : 'en';
      const tips = {
        en: [
          '💡 Start with a preset to see a working example!',
          '💡 Hover over field labels to see helpful tooltips',
          '💡 Lower doses = fewer side effects, but slower tumor reduction',
          '💡 Check the "Healthy Cells" meter - keep it green!',
          '💡 Drug interactions matter - watch the interaction strength',
          '💡 Try the Help button (?) to learn about each parameter',
          '💡 Compare different time horizons to see long-term effects',
        ],
        it: [
          '💡 Inizia con un preset per vedere un esempio funzionante!',
          '💡 Passa sopra le etichette dei campi per vedere suggerimenti utili',
          '💡 Dosi più basse = meno effetti collaterali, ma riduzione tumorale più lenta',
          '💡 Controlla il metro "Cellule Sane" - mantienilo verde!',
          '💡 Le interazioni tra farmaci contano - guarda la forza dell\'interazione',
          '💡 Prova il pulsante Aiuto (?) per saperne di più su ogni parametro',
          '💡 Confronta diversi orizzonti temporali per vedere gli effetti a lungo termine',
        ],
      };

      const langTips = tips[lang] || tips.en;
      return langTips[Math.floor(Math.random() * langTips.length)];
    }

    reset() {
      this.data = this.getDefaultData();
      this.save();
    }
  }

  // 🎨 UI Helper Functions
  function showAchievement(achievement) {
    if (!achievement) return;
    
    const lang = window.getCurrentLang ? window.getCurrentLang() : 'en';
    const title = achievement.title[lang] || achievement.title.en;
    const description = achievement.description[lang] || achievement.description.en;
    
    const html = `
      <div class="achievement-toast">
        <div class="achievement-icon">${achievement.icon}</div>
        <div class="achievement-content">
          <div class="achievement-title">🏆 ${title}</div>
          <div class="achievement-description">${description}</div>
          <div class="achievement-points">+${achievement.points} points</div>
        </div>
      </div>
    `;
    
    const container = document.createElement('div');
    container.innerHTML = html;
    container.className = 'achievement-notification';
    document.body.appendChild(container);
    
    setTimeout(function () {
      container.classList.add('show');
    }, 100);
    
    setTimeout(function () {
      container.classList.remove('show');
      setTimeout(function () {
        document.body.removeChild(container);
      }, 500);
    }, 5000);
  }

  function showProgressBadge() {
    const tracker = window.MMPortal.progress;
    if (!tracker) return;
    
    const level = tracker.getLevel();
    const progress = tracker.getProgress();
    const lang = window.getCurrentLang ? window.getCurrentLang() : 'en';
    
    const badge = document.getElementById('progress-badge');
    if (badge) {
      badge.innerHTML = `
        <div class="level-badge">
          <span class="level-icon">${level.icon}</span>
          <span class="level-name">${level.name[lang]}</span>
          <span class="level-points">${tracker.data.points} pts</span>
        </div>
        <div class="progress-bar-container">
          <div class="progress-bar-fill" style="width: ${progress.percentage}%"></div>
        </div>
      `;
    }
  }

  // 🚀 Initialize
  function init() {
    if (!root.MMPortal) {
      root.MMPortal = {};
    }
    root.MMPortal.progress = new ProgressTracker();
    root.showAchievement = showAchievement;
    root.showProgressBadge = showProgressBadge;
    
    // Show progress badge on load
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', showProgressBadge);
    } else {
      showProgressBadge();
    }
  }

  init();
})(typeof window !== 'undefined' ? window : globalThis);
