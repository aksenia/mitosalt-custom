# SAltShaker Docker image
# Minimal container for mitochondrial structural variant visualization
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install only Python and essential tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python-is-python3 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/saltshaker

# Copy saltshaker package source code
COPY saltshaker/ /opt/saltshaker/saltshaker/
COPY setup.py pyproject.toml README.md /opt/saltshaker/

# Install saltshaker package and its Python dependencies
# Dependencies from setup.py: pandas, numpy, matplotlib, biopython
RUN pip install --no-cache-dir --break-system-packages .

# Set working directory for user data
WORKDIR /data

# Add saltshaker to PATH (already in PATH via pip install)
ENV PYTHONPATH="/opt/saltshaker:${PYTHONPATH}"

# Set entrypoint to saltshaker command
ENTRYPOINT ["saltshaker"]

# Default command (can be overridden)
CMD ["--help"]