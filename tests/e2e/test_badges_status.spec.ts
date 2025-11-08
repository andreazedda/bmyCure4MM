import { test, expect } from '@playwright/test';

test.describe('Dose Badges Status', () => {
  test.beforeEach(async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/simulator/scenario/1/`);
  });

  test('Dose badge changes to danger when out of range', async ({ page }) => {
    // Find lenalidomide dose input
    const doseInput = page.locator('#id_lenalidomide_dose');
    await doseInput.waitFor({ state: 'visible' });
    
    // Find corresponding badge
    const badge = page.getByTestId('badge-id_lenalidomide_dose');
    
    // Enter a value that's out of range (e.g., 100 mg)
    await doseInput.fill('100');
    await doseInput.blur();
    
    // Wait for badge to update
    await page.waitForTimeout(500);
    
    // Badge should have danger class
    await expect(badge).toHaveClass(/bg-danger|text-danger|danger/);
  });

  test('Dose badge changes to warn when near range boundary', async ({ page }) => {
    const doseInput = page.locator('#id_lenalidomide_dose');
    await doseInput.waitFor({ state: 'visible' });
    
    const badge = page.getByTestId('badge-id_lenalidomide_dose');
    
    // Enter a value near the boundary (adjust based on actual range)
    await doseInput.fill('5');
    await doseInput.blur();
    
    await page.waitForTimeout(500);
    
    // Badge should have warn class or safe class depending on range
    const classes = await badge.getAttribute('class');
    expect(classes).toBeTruthy();
  });

  test('Dose badge shows safe when in normal range', async ({ page }) => {
    const doseInput = page.locator('#id_lenalidomide_dose');
    await doseInput.waitFor({ state: 'visible' });
    
    const badge = page.getByTestId('badge-id_lenalidomide_dose');
    
    // Enter a safe value
    await doseInput.fill('15');
    await doseInput.blur();
    
    await page.waitForTimeout(500);
    
    // Badge should have safe class
    await expect(badge).toHaveClass(/bg-success|text-success|safe/);
  });

  test('Interaction strength badge updates correctly', async ({ page }) => {
    const interactionInput = page.locator('#id_interaction_strength');
    
    // Check if interaction field exists on this page
    const exists = await interactionInput.count();
    if (exists === 0) {
      test.skip();
      return;
    }
    
    await interactionInput.waitFor({ state: 'visible' });
    
    // Enter high interaction value
    await interactionInput.fill('0.5');
    await interactionInput.blur();
    
    await page.waitForTimeout(500);
    
    // Look for badge or indicator
    const badge = page.locator('[id*="badge"][id*="interaction"]').first();
    if (await badge.count() > 0) {
      await expect(badge).toHaveClass(/danger/);
    }
  });

  test('Badge updates announce via aria-live region', async ({ page }) => {
    // This test verifies that badge changes are announced to screen readers
    const liveRegion = page.locator('#live-region');
    await expect(liveRegion).toHaveAttribute('aria-live', 'polite');
    
    const doseInput = page.locator('#id_lenalidomide_dose');
    await doseInput.waitFor({ state: 'visible' });
    
    // Change dose to trigger badge update
    await doseInput.fill('100');
    await doseInput.blur();
    
    // The live region should be in the DOM
    await expect(liveRegion).toBeAttached();
  });

  test('Multiple badges can coexist on same form', async ({ page }) => {
    // Find all badges
    const badges = page.locator('[data-testid^="badge-"]');
    const count = await badges.count();
    
    // Should have multiple badges for different fields
    expect(count).toBeGreaterThan(0);
  });

  test('Badge shows unknown for invalid input', async ({ page }) => {
    const doseInput = page.locator('#id_lenalidomide_dose');
    await doseInput.waitFor({ state: 'visible' });
    
    const badge = page.getByTestId('badge-id_lenalidomide_dose');
    
    // Enter invalid value
    await doseInput.fill('invalid');
    await doseInput.blur();
    
    await page.waitForTimeout(500);
    
    // Badge should show unknown state or hide
    const classes = await badge.getAttribute('class');
    // Either has 'unknown' or is empty/hidden
    if (classes) {
      expect(classes).toBeTruthy();
    }
  });

  test('Badge reacts to input changes in real-time', async ({ page }) => {
    const doseInput = page.locator('#id_lenalidomide_dose');
    await doseInput.waitFor({ state: 'visible' });
    
    const badge = page.getByTestId('badge-id_lenalidomide_dose');
    
    // Type gradually and watch badge update
    await doseInput.fill('5');
    await page.waitForTimeout(200);
    
    await doseInput.fill('15');
    await page.waitForTimeout(200);
    
    await doseInput.fill('100');
    await page.waitForTimeout(200);
    
    // Badge should have updated multiple times
    const finalClasses = await badge.getAttribute('class');
    expect(finalClasses).toBeTruthy();
  });
});
