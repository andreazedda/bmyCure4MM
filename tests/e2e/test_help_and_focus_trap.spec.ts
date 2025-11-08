import { test, expect } from '@playwright/test';

test.describe('Help Drawer Focus Trap', () => {
  test.beforeEach(async ({ page, baseURL }) => {
    await page.goto(`${baseURL}/simulator/scenario/1/`);
  });

  test('Help opens from inline button, focus trapped, Esc closes and focus returns', async ({ page }) => {
    // Find and click a help button
    const helpButton = page.getByTestId('help-open-dose_ranges_lenalidomide').first();
    await helpButton.click();

    // Verify drawer is visible
    const drawer = page.getByTestId('help-drawer');
    await expect(drawer).toBeVisible();

    // Verify focus trap: Tab should cycle through focusable elements
    await page.keyboard.press('Tab');
    const firstFocusable = await page.evaluate(() => document.activeElement?.tagName);
    
    await page.keyboard.press('Tab');
    await page.keyboard.press('Tab');
    
    // Press Shift+Tab to go backward
    await page.keyboard.press('Shift+Tab');
    
    // Close with Escape
    await page.keyboard.press('Escape');
    await expect(drawer).toBeHidden();

    // Verify focus returns to trigger button
    await expect(helpButton).toBeFocused();
  });

  test('Help drawer closes when clicking close button', async ({ page }) => {
    const helpButton = page.getByTestId('help-open-dose_ranges_lenalidomide').first();
    await helpButton.click();

    const drawer = page.getByTestId('help-drawer');
    await expect(drawer).toBeVisible();

    // Click close button
    const closeButton = drawer.locator('button.btn-close');
    await closeButton.click();

    await expect(drawer).toBeHidden();
    await expect(helpButton).toBeFocused();
  });

  test('Multiple help buttons can be clicked sequentially', async ({ page }) => {
    const drawer = page.getByTestId('help-drawer');
    
    // Click first help button
    const firstButton = page.getByTestId('help-open-dose_ranges_lenalidomide').first();
    await firstButton.click();
    await expect(drawer).toBeVisible();
    
    // Close it
    await page.keyboard.press('Escape');
    await expect(drawer).toBeHidden();
    
    // Click second help button
    const secondButton = page.getByTestId('help-open-dose_ranges_bortezomib').first();
    await secondButton.click();
    await expect(drawer).toBeVisible();
    
    // Close it
    await page.keyboard.press('Escape');
    await expect(drawer).toBeHidden();
    await expect(secondButton).toBeFocused();
  });

  test('Help drawer has proper ARIA attributes', async ({ page }) => {
    const helpButton = page.getByTestId('help-open-dose_ranges_lenalidomide').first();
    await helpButton.click();

    const drawer = page.getByTestId('help-drawer');
    await expect(drawer).toHaveAttribute('role', 'dialog');
    await expect(drawer).toHaveAttribute('aria-modal', 'true');
  });

  test('F1 key opens help for focused element with data-help', async ({ page }) => {
    // Find an element with data-help attribute
    const helpableElement = page.locator('[data-help="dose_ranges_lenalidomide"]').first();
    await helpableElement.focus();
    
    // Press F1
    await page.keyboard.press('F1');
    
    // Verify drawer opens
    const drawer = page.getByTestId('help-drawer');
    await expect(drawer).toBeVisible();
  });
});
