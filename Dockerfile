FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy project files
COPY . .

# Create necessary directories
RUN mkdir -p logs media static

# Collect static files (for production)
RUN python manage.py collectstatic --noinput

EXPOSE 8001

ENTRYPOINT ["sh", "/app/entrypoint.sh"]

# Default command (dev). Production should override to gunicorn.
CMD ["python", "manage.py", "runserver", "0.0.0.0:8001"]
