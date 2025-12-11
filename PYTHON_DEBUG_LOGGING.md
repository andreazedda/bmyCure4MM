# 🐍 Python Debug Logging - Implementation Summary

## Overview
Added **comprehensive scrupulous debug logging** to all Python components involved in the binding visualizer workflow: Django views, API client, and the binding_visualizer.py script itself.

---

## 📦 Files Enhanced with Logging

### 1. **chemtools/views.py** (Django Views)
Tracks HTTP requests, form processing, job creation, and metadata enrichment.

**Functions Enhanced:**
- `binding_viz()` - Form submission and job creation
- `job_detail()` - Results display and API enrichment

**What's Logged:**
```python
[DEBUG] 🌐 binding_viz view called - Method: POST, User: admin
[DEBUG] 📝 Form validated - PDB ID: 5LF3, Ligand: BTZ
[DEBUG] 📊 API Preferences captured: {'fetch_validation': True, ...}
[DEBUG] ✅ Enabled APIs (7): fetch_validation, fetch_interactions, ...
[DEBUG] 💾 Job created - ID: 42, Kind: BIND, Status: Q
[DEBUG] 🔗 Job redirect URL: /chemtools/job/42/
[DEBUG] 🚀 Attempting to enqueue job 42...
[DEBUG] ✅ Job 42 queued successfully (Celery)
```

### 2. **chemtools/pdb_api_client.py** (API Client)
Tracks all external API calls with timing, status codes, and response sizes.

**Methods Enhanced:**
- `get_pdb_summary()`
- `get_uniprot_info()`
- (All other API methods follow same pattern)

**What's Logged:**
```python
[DEBUG] 🌐 API CALL: get_pdb_summary('5LF3')
[DEBUG]   URL: https://data.rcsb.org/rest/v1/core/entry/5LF3
[DEBUG]   Timeout: 10s
[DEBUG]   Status: 200
[DEBUG]   Response time: 0.34s
[DEBUG]   Response size: 15234 bytes
[DEBUG] ✅ Successfully fetched PDB summary for 5LF3 in 0.34s
[DEBUG]   Data keys: ['entry', 'rcsb_primary_citation', 'struct', ...]
```

### 3. **modules/binding_visualizer/binding_visualizer.py** (Core Script)
Tracks script initialization, PDB data fetching, and structure processing.

**What's Logged:**
```python
══════════════════════════════════════════════════════════════════
  bmyCure4MM - BINDING VISUALIZER - DEBUG MODE ACTIVE
══════════════════════════════════════════════════════════════════
🔍 [DEBUG] All operations will be logged in detail
📊 [DEBUG] API calls and data processing tracked
⏱️  [DEBUG] Performance timings enabled
🔬 [DEBUG] Data extraction and validation logged
══════════════════════════════════════════════════════════════════

[DEBUG] 🚀 Starting script...
[DEBUG] 📅 Timestamp: 2025-12-11T10:30:45.123456
[DEBUG] 🐍 Python version: 3.9.6 (default, ...)
[DEBUG] 💻 Platform: macOS-14.0-arm64
[DEBUG] 👤 User: andrea
[DEBUG] 📁 Script folder: /Volumes/nvme/Github/bmyCure4MM/modules/binding_visualizer
[DEBUG] 📝 Loading configuration from .../binding_visualizer.yaml
[DEBUG] ✅ Config file exists (size: 2345 bytes)

══════════════════════════════════════════════════════════════════
  STARTING MAIN WORKFLOW
══════════════════════════════════════════════════════════════════
[DEBUG] 🎯 Entered main()
[DEBUG] 📋 Log file: .../binding_visualizer.log
[DEBUG] 📊 Configuration keys: ['pdb_id', 'viewer', 'visualization', ...]
[DEBUG] 🔬 PDB ID: 5LF3
[DEBUG] 📐 Viewer dimensions: 800x600

[DEBUG] 🌐 Fetching PDB data for 5LF3...
[DEBUG] ✅ PDB data fetched in 0.45s
[DEBUG] 📄 PDB data size: 245678 bytes
[DEBUG] 📄 PDB data lines: 3456 lines
```

---

## 📊 Log Output Locations

### Django Logs
Django logs go to the console (terminal where `python manage.py runserver` runs) and can be configured in `settings.py`.

**View Django logs:**
```bash
# In the terminal running the dev server
python manage.py runserver
# Logs appear here automatically
```

### Binding Visualizer Logs
The binding_visualizer.py script creates its own log file:

**Location:** `modules/binding_visualizer/binding_visualizer.log`

**View the log:**
```bash
# Tail the log file to see real-time updates
tail -f /Volumes/nvme/Github/bmyCure4MM/modules/binding_visualizer/binding_visualizer.log

# Or view the entire file
cat /Volumes/nvme/Github/bmyCure4MM/modules/binding_visualizer/binding_visualizer.log
```

### API Client Logs
API client logs appear in Django's logger output (console).

---

## 🔍 What Each Component Tracks

### Django Views (`views.py`)

| Event | Log Level | What's Tracked |
|-------|-----------|----------------|
| View entry | INFO | Method, user, timestamp |
| Form validation | INFO | PDB ID, ligand, validation result |
| API preferences | DEBUG | All 10 preference flags |
| Job creation | INFO | Job ID, kind, status |
| Job enqueueing | INFO | Queued/failed status |
| HTML loading | DEBUG/INFO | File size, line count |
| CSV loading | DEBUG/INFO | Row count, columns |
| Metadata enrichment | INFO | PDB ID, timing |
| View completion | INFO | Total elapsed time |

### API Client (`pdb_api_client.py`)

| Event | Log Level | What's Tracked |
|-------|-----------|----------------|
| API call start | INFO | Method, PDB ID/UniProt ID |
| Request details | DEBUG | URL, timeout |
| Response received | DEBUG | Status code, response time, size |
| Success | INFO | Total time, data keys |
| Failure | WARNING | Error type, elapsed time |

### Binding Visualizer (`binding_visualizer.py`)

| Event | Log Level | What's Tracked |
|-------|-----------|----------------|
| Script start | INFO | Timestamp, Python version, platform, user |
| Config loading | DEBUG | File path, existence, size |
| Main workflow | INFO | Entry point, config keys |
| PDB fetch | INFO | PDB ID, fetch time, data size |
| Parsing | INFO | Metadata extraction |
| Visualization | INFO | Output file creation |
| PDF generation | INFO | Success/failure, file path |
| Errors | ERROR | Full traceback |

---

## 🎯 How to Use Python Debugging

### 1. Start Django Development Server with Logging
```bash
cd /Volumes/nvme/Github/bmyCure4MM
source venv/bin/activate
python manage.py runserver

# Logs will appear in this terminal
```

### 2. Submit a Binding Visualizer Job
1. Navigate to: http://127.0.0.1:8000/chemtools/binding-viz/
2. Select an MM drug preset (e.g., "🔴 Bortezomib")
3. Check API sources to enable
4. Submit the form

**Watch the terminal** - you'll see:
```
[DEBUG] 🌐 binding_viz view called - Method: POST, User: admin
[DEBUG] 📝 Form validated - PDB ID: 5LF3, Ligand: BTZ
[DEBUG] 📊 API Preferences captured: {...}
[DEBUG] ✅ Enabled APIs (7): fetch_validation, fetch_interactions, ...
[DEBUG] 💾 Job created - ID: 42, Kind: BIND, Status: Q
[DEBUG] 🚀 Attempting to enqueue job 42...
[DEBUG] ✅ Job 42 queued successfully (Celery)
```

### 3. View Job Detail Page
Navigate to the job detail page (e.g., `/chemtools/job/42/`)

**Watch the terminal** - you'll see:
```
[DEBUG] 📋 job_detail view called - Job ID: 42, User: admin
[DEBUG] 💼 Job loaded - Kind: BIND, Status: C, Created: 2025-12-11 10:30:45
[DEBUG]   Input A: 5LF3, Input B: BTZ
[DEBUG]   Has HTML: True, Has CSV: False
[DEBUG] 📄 Loading HTML output...
[DEBUG] ✅ HTML loaded - Size: 372009 chars, Lines: 1234
[DEBUG] 🔬 Processing binding job metadata for PDB: 5LF3
[DEBUG]   Basic URLs configured: ['pdb_id', 'rcsb_url', ...]
[DEBUG] 🔧 API preferences loaded: 7 sources enabled
[DEBUG] 🌐 Enriching metadata with API data...
[DEBUG] 🌐 API CALL: get_pdb_summary('5LF3')
[DEBUG]   URL: https://data.rcsb.org/rest/v1/core/entry/5LF3
[DEBUG]   Status: 200
[DEBUG]   Response time: 0.34s
[DEBUG] ✅ Successfully fetched PDB summary for 5LF3 in 0.34s
[DEBUG] ✅ Metadata enrichment completed in 2.45s
[DEBUG] ⏱️ job_detail view completed in 2.67s
```

### 4. Check Binding Visualizer Log File
```bash
# View real-time logs
tail -f /Volumes/nvme/Github/bmyCure4MM/modules/binding_visualizer/binding_visualizer.log

# Search for specific PDB ID
grep "5LF3" /Volumes/nvme/Github/bmyCure4MM/modules/binding_visualizer/binding_visualizer.log

# View last 100 lines
tail -n 100 /Volumes/nvme/Github/bmyCure4MM/modules/binding_visualizer/binding_visualizer.log
```

---

## 🐛 Troubleshooting with Python Logs

### Issue: Job Stuck in Queue
**Check:**
```
# Look for job creation
grep "Job created - ID: 42" terminal_output

# Check if job was enqueued
grep "Job 42 queued" terminal_output

# If not queued, check for errors
grep "ERROR" terminal_output
```

### Issue: API Not Returning Data
**Check:**
```
# Look for API calls
grep "API CALL:" terminal_output

# Check response status
grep "Status:" terminal_output

# Look for failures
grep "Failed to fetch" terminal_output
```

### Issue: Slow Performance
**Check timing logs:**
```
# View all timing information
grep "in [0-9.]*s" terminal_output

# Compare:
# - PDB data fetched in: X.XXs
# - Metadata enrichment completed in: X.XXs
# - job_detail view completed in: X.XXs
```

### Issue: Missing Data
**Check:**
```
# See what was loaded
grep "loaded -" terminal_output

# Check enrichment
grep "Enriched data keys:" terminal_output

# Look for warnings
grep "⚠️" terminal_output
```

---

## 📈 Performance Benchmarks

**Expected timings:**
- Django view entry to response: 0.1-5s
- PDB data fetch: 0.3-2s  
- Single API call: 0.2-1s
- Metadata enrichment (all APIs): 2-10s
- HTML file loading: <0.1s

**If significantly slower:**
- Check network connectivity
- Check API rate limits
- Check file system performance
- Look for timeout errors in logs

---

## 🔧 Customizing Log Verbosity

### Change Django Logging Level
Edit `mmportal/settings.py`:

```python
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
        },
    },
    'loggers': {
        'chemtools': {
            'handlers': ['console'],
            'level': 'DEBUG',  # Change to INFO, WARNING, or ERROR
        },
    },
}
```

### Change Binding Visualizer Logging
In `binding_visualizer.py`, change:
```python
logging.basicConfig(
    level=logging.DEBUG,  # Change to INFO, WARNING, or ERROR
    ...
)
```

---

## 📚 Log Level Guide

| Level | When to Use | What It Shows |
|-------|-------------|---------------|
| **DEBUG** | Development, troubleshooting | Everything (very verbose) |
| **INFO** | Normal operation | Key events and milestones |
| **WARNING** | Potential issues | Non-fatal problems |
| **ERROR** | Failures | Fatal errors and exceptions |

---

## ✅ Verification Checklist

Test Python logging is working:

1. ✅ Start dev server, see Django startup logs
2. ✅ Submit binding viz form, see view logs
3. ✅ Check terminal for job creation logs
4. ✅ Navigate to job detail, see loading logs
5. ✅ See API call logs with timings
6. ✅ Check binding_visualizer.log exists
7. ✅ Tail log file, see real-time updates

---

**Complete logging coverage achieved!** 🎉

Python logs track:
- ✅ HTTP requests and responses
- ✅ Form validation and processing
- ✅ Database operations
- ✅ All external API calls
- ✅ File I/O operations
- ✅ Performance timings
- ✅ Error conditions with full context
