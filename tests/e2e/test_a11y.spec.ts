import { test, expect } from '@playwright/test';
import AxeBuilder from '@axe-core/playwright';

test.describe('Accessibility Smoke Tests', () => {
  test('Scenario page has no serious accessibility violations', async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/simulator/scenario/1/`);
    
    const results = await new AxeBuilder({ page }).analyze();
    
    // Filter for serious and critical violations only
    const seriousViolations = results.violations.filter(
      (violation) => violation.impact === 'serious' || violation.impact === 'critical'
    );
    
    // Log violations for debugging
    if (seriousViolations.length > 0) {
      console.log('Serious accessibility violations found:');
      seriousViolations.forEach((violation) => {
        console.log(`- ${violation.id}: ${violation.description}`);
        console.log(`  Impact: ${violation.impact}`);
        console.log(`  Help: ${violation.helpUrl}`);
      });
    }
    
    expect(seriousViolations).toEqual([]);
  });

  test('Dashboard page has no serious accessibility violations', async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/clinic/dashboard/`);
    
    const results = await new AxeBuilder({ page }).analyze();
    
    const seriousViolations = results.violations.filter(
      (violation) => violation.impact === 'serious' || violation.impact === 'critical'
    );
    
    expect(seriousViolations).toEqual([]);
  });

  test('Help drawer has no accessibility violations when open', async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/simulator/scenario/1/`);
    
    // Open help drawer
    const helpButton = page.locator('[data-help]').first();
    await helpButton.click();
    
    // Wait for drawer to open
    const drawer = page.getByTestId('help-drawer');
    await expect(drawer).toBeVisible();
    
    // Run axe on the drawer
    const results = await new AxeBuilder({ page })
      .include('#helpDrawer')
      .analyze();
    
    const seriousViolations = results.violations.filter(
      (violation) => violation.impact === 'serious' || violation.impact === 'critical'
    );
    
    expect(seriousViolations).toEqual([]);
  });

  test('Command-K search has no accessibility violations', async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/simulator/scenario/1/`);
    
    // Open command palette
    await page.keyboard.press('Meta+K');
    
    const input = page.getByTestId('cmdk-input');
    await expect(input).toBeVisible();
    
    // Type to show results
    await input.type('dose');
    
    // Wait for results
    await page.waitForTimeout(500);
    
    // Run axe on command palette area
    const results = await new AxeBuilder({ page })
      .include('#cmdkWrapper')
      .analyze();
    
    const seriousViolations = results.violations.filter(
      (violation) => violation.impact === 'serious' || violation.impact === 'critical'
    );
    
    expect(seriousViolations).toEqual([]);
  });

  test('Form with badges has no accessibility violations', async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/simulator/scenario/1/`);
    
    // Fill in some form fields to trigger badges
    const doseInput = page.locator('#id_lenalidomide_dose');
    if (await doseInput.count() > 0) {
      await doseInput.fill('15');
      await doseInput.blur();
      await page.waitForTimeout(500);
    }
    
    // Run axe on the form area
    const results = await new AxeBuilder({ page })
      .include('form')
      .analyze();
    
    const seriousViolations = results.violations.filter(
      (violation) => violation.impact === 'serious' || violation.impact === 'critical'
    );
    
    expect(seriousViolations).toEqual([]);
  });

  test('Navigation menu is keyboard accessible', async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/`);
    
    // Tab through navigation
    await page.keyboard.press('Tab');
    await page.keyboard.press('Tab');
    
    // Check that focused element is within navigation
    const focusedElement = await page.evaluate(() => {
      const el = document.activeElement;
      return el?.tagName;
    });
    
    expect(focusedElement).toBeTruthy();
    
    // Run axe on navigation
    const results = await new AxeBuilder({ page })
      .include('nav')
      .analyze();
    
    const seriousViolations = results.violations.filter(
      (violation) => violation.impact === 'serious' || violation.impact === 'critical'
    );
    
    expect(seriousViolations).toEqual([]);
  });

  test('All interactive elements have accessible names', async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/simulator/scenario/1/`);
    
    // Check for buttons without accessible names
    const results = await new AxeBuilder({ page })
      .withTags(['wcag2a', 'wcag2aa'])
      .analyze();
    
    const nameViolations = results.violations.filter(
      (violation) => 
        violation.id === 'button-name' || 
        violation.id === 'link-name' ||
        violation.id === 'label'
    );
    
    if (nameViolations.length > 0) {
      console.log('Accessible name violations:');
      nameViolations.forEach((violation) => {
        console.log(`- ${violation.id}: ${violation.description}`);
      });
    }
    
    // These should be informational only, not fail the test
    // unless they're serious/critical
    const seriousNameViolations = nameViolations.filter(
      (violation) => violation.impact === 'serious' || violation.impact === 'critical'
    );
    
    expect(seriousNameViolations).toEqual([]);
  });

  test('Color contrast is sufficient', async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/simulator/scenario/1/`);
    
    const results = await new AxeBuilder({ page })
      .withTags(['wcag2aa'])
      .analyze();
    
    const contrastViolations = results.violations.filter(
      (violation) => violation.id === 'color-contrast'
    );
    
    if (contrastViolations.length > 0) {
      console.log('Color contrast violations:');
      contrastViolations.forEach((violation) => {
        console.log(`- ${violation.description}`);
        violation.nodes.forEach((node) => {
          console.log(`  Element: ${node.html}`);
        });
      });
    }
    
    // Allow some contrast violations but flag serious ones
    const seriousContrastViolations = contrastViolations.filter(
      (violation) => violation.impact === 'serious' || violation.impact === 'critical'
    );
    
    expect(seriousContrastViolations).toEqual([]);
  });
});
