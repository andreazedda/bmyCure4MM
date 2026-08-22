FROM ghcr.io/astral-sh/uv:0.12.3 AS uv
FROM python:3.11.11-slim

WORKDIR /app

ENV PATH="/app/.venv/bin:$PATH" \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

COPY --from=uv /uv /uvx /bin/

# Install the exact application graph before copying source for layer reuse.
COPY pyproject.toml uv.lock .python-version ./
RUN uv sync --frozen --no-dev --extra chemistry --no-install-project

# Copy project files
COPY . .

# Create necessary directories
RUN mkdir -p logs media static

# Collect static files without baking a build-time log or fallback secret into
# the image. These values exist only for this deterministic build command.
RUN DJANGO_DEBUG=0 \
    DJANGO_SECRET_KEY=container-build-only-synthetic-secret \
    DJANGO_LOGS_ROOT=/tmp/bmycure4mm-build-logs \
    python manage.py collectstatic --noinput \
    && rm -rf /tmp/bmycure4mm-build-logs

EXPOSE 8001

ENTRYPOINT ["sh", "/app/entrypoint.sh"]

# Default command (dev). Production should override to gunicorn.
CMD ["python", "manage.py", "runserver", "0.0.0.0:8001"]
