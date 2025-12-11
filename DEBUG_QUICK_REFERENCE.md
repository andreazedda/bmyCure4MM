# 🔍 Quick Debug Reference

## Open Console
**Mac**: `Cmd + Option + J` (Chrome) or `Cmd + Option + C` (Safari)
**Windows**: `F12` or `Ctrl + Shift + J`

## What to Look For

### ✅ Success Indicators
```
✅ window.viewer already available
✅ Content reorganization completed
✅ Controls reorganization completed
✅ Enabled 16 of 16 control buttons
🏁 Initialization complete: { success: true }
📚 All debugging enabled. Click any button to see detailed logs.
```

### ❌ Problem Indicators
```
❌ No container found
❌ No info div found
❌ window.viewer not found
⚠️ Button is disabled, click blocked!
⏳ Viewer not ready, will retry...
```

## Console Filters

Type these in the console filter box to focus on specific areas:

| Filter | Shows |
|--------|-------|
| `[DEBUG]` | All debug logs |
| `Button Click` | Button interactions |
| `Viewer` | Viewer state checks |
| `Extracting` | Data extraction |
| `ms` | Performance timings |
| `enabled` | Button enable operations |
| `❌` or `⚠️` | Errors/warnings only |
| `✅` | Success messages only |

## Quick Tests

### Test 1: Check Initialization
1. Open console (F12)
2. Refresh page
3. Look for: `🏁 Initialization complete: { success: true }`

### Test 2: Check Buttons
1. Scroll to "Final State Summary"
2. Verify: `enabled: 16` and `disabled: 0`

### Test 3: Test Click
1. Click "Cartoon" button
2. Look for: `[DEBUG] 🖱️ Button Click #0`
3. Verify: `✅ Button enabled, onclick should fire`

### Test 4: Check Viewer
1. Expand "Viewer State" object
2. Verify all properties are `true`

## Common Issues & Solutions

### Issue: Buttons Disabled
**Console shows**: `disabled: 16`
**Solution**: Check if initialization completed, may need to wait

### Issue: Viewer Not Found
**Console shows**: `❌ window.viewer not found`
**Solution**: py3Dmol library may not have loaded, check Network tab

### Issue: No Click Response
**Console shows**: `Button Click` but no style change
**Solution**: Check if `window.viewer.setStyle` exists in Viewer State

### Issue: Missing Data
**Console shows**: `extractField("X"): (not found)`
**Solution**: binding_visualizer.py may not have generated that field

## Debug Data Export

To save console output for analysis:
1. Right-click in console
2. Select "Save as..."
3. Save as .log file
4. Share with developer

## Performance Benchmarks

**Expected values**:
- Initialization: 100-1000ms
- Content reorganization: 50-300ms  
- Button click response: <10ms
- Viewer render: 100-500ms

If times are significantly higher, check for:
- Large molecule structures (>10k atoms)
- Slow network (API data fetching)
- Browser performance issues

---

**Last Updated**: December 11, 2025
