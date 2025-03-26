# Use the official Rocker R base image
FROM rocker/r-base:latest

# Metadata
LABEL maintainer="your.email@example.com"
LABEL description="Docker container for development of InSilicoVA-sim (R, Java, and required packages)"

# Install OpenJDK 11 and system dependencies
RUN apt-get update && apt-get install -y \
    openjdk-11-jdk \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set JAVA_HOME and LD_LIBRARY_PATH environment variables
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
ENV LD_LIBRARY_PATH=$JAVA_HOME/lib/server

# Reconfigure R to use the new Java settings
RUN R CMD javareconf

# Install required R packages (rJava, openVA, MCMCpack, xtable)
RUN R -e "install.packages(c('rJava', 'openVA', 'MCMCpack', 'xtable'), repos='https://cran.rstudio.com')"
COPY . /InSilicoVA-sim

# Declare a volume so that you can mount your working directory for development
VOLUME ["/InSilicoVA-sim"]

# Set the working directory
WORKDIR /InSilicoVA-sim

# Default command: start an interactive bash shell
CMD ["bash"]