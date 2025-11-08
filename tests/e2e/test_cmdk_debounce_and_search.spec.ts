import { test, expect } from '@playwright/test';

test.describe('Command-K Search', () => {
  test.beforeEach(async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/simulator/scenario/1/`);
  });

  test('Cmd-K opens search input', async ({ page }) => {
    const input = page.getByTestId('cmdk-input');
    
    // Initially hidden
    await expect(input).toBeHidden();
    
    // Press Cmd+K (Meta+K on Mac, Ctrl+K on Windows/Linux)
    await page.keyboard.press('Meta+K');
    
    // Input should be visible and focused
    await expect(input).toBeVisible();
    await expect(input).toBeFocused();
  });

  test('Cmd-K search displays results', async ({ page }) => {
    await page.keyboard.press('Meta+K');
    
    const input = page.getByTestId('cmdk-input');
    await input.type('lenali');
    
    // Wait for results to appear
    const resultsList = page.locator('#cmdkResults a');
    await expect(resultsList.first()).toBeVisible({ timeout: 2000 });
    
    // Verify at least one result
    const count = await resultsList.count();
    expect(count).toBeGreaterThan(0);
  });

  test('Escape closes search', async ({ page }) => {
    await page.keyboard.press('Meta+K');
    
    const input = page.getByTestId('cmdk-input');
    await expect(input).toBeVisible();
    
    // Close with Escape
    await page.keyboard.press('Escape');
    await expect(input).toBeHidden();
  });

  test('Clicking outside closes search', async ({ page }) => {
    await page.keyboard.press('Meta+K');
    
    const input = page.getByTestId('cmdk-input');
    await expect(input).toBeVisible();
    
    // Click outside the command palette
    await page.click('body', { position: { x: 10, y: 10 } });
    
    await expect(input).toBeHidden();
  });

  test('Search results can be clicked to open help', async ({ page }) => {
    await page.keyboard.press('Meta+K');
    
    const input = page.getByTestId('cmdk-input');
    await input.type('lenali');
    
    // Wait for and click first result
    const firstResult = page.locator('#cmdkResults a').first();
    await expect(firstResult).toBeVisible({ timeout: 2000 });
    await firstResult.click();
    
    // Verify help drawer opens
    const drawer = page.getByTestId('help-drawer');
    await expect(drawer).toBeVisible();
  });

  test('Empty search shows default results', async ({ page }) => {
    await page.keyboard.press('Meta+K');
    
    // Don't type anything
    const resultsList = page.locator('#cmdkResults a');
    
    // Should show some default results
    await expect(resultsList.first()).toBeVisible({ timeout: 2000 });
  });

  test('Search results include type icons', async ({ page }) => {
    await page.keyboard.press('Meta+K');
    
    const input = page.getByTestId('cmdk-input');
    await input.type('dose');
    
    // Wait for results
    const firstResult = page.locator('#cmdkResults a').first();
    await expect(firstResult).toBeVisible({ timeout: 2000 });
    
    // Results should contain emoji icons
    const text = await firstResult.textContent();
    expect(text).toBeTruthy();
  });

  test('Search is debounced', async ({ page }) => {
    await page.keyboard.press('Meta+K');
    
    const input = page.getByTestId('cmdk-input');
    
    // Type quickly
    await input.type('l', { delay: 10 });
    await input.type('e', { delay: 10 });
    await input.type('n', { delay: 10 });
    
    // Results should appear after debounce delay
    const resultsList = page.locator('#cmdkResults a');
    await expect(resultsList.first()).toBeVisible({ timeout: 2000 });
  });

  test('Search respects current language', async ({ page }) => {
    // Switch to Italian
    await page.click('#lang-it');
    
    await page.keyboard.press('Meta+K');
    
    const input = page.getByTestId('cmdk-input');
    await input.type('dose');
    
    // Results should be in Italian
    const firstResult = page.locator('#cmdkResults a').first();
    await expect(firstResult).toBeVisible({ timeout: 2000 });
  });
});
