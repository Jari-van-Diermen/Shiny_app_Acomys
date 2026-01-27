# syntax=docker/dockerfile:1
FROM rocker/shiny:4.5.2

### Install R packages
RUN R -e "install.packages(c('renv', 'devtools'), repos = c(CRAN = 'https://cloud.r-project.org'))"

# Copy renv files to container
WORKDIR /project
COPY renv.lock renv.lock

# Set up renv environment and indicate library path
RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Install required packages
RUN apt-get update -qq && \
    apt-get install -y openssl
RUN apt-get update -qq && \
    apt-get install -y libssl-dev \
                       libxml2-dev \
                       libharfbuzz-dev \
                       libfribidi-dev \
                       libfreetype6-dev \
                       libfreetype-dev \
                       libpng-dev \
                       libtiff5-dev \
                       libjpeg-dev

# restore renv environment (i.e. install R packages)
RUN R -e "renv::restore()"

### Retrieve Shiny app files
COPY Shiny_Acomys Shiny_Acomys

### Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Expose the application port
EXPOSE 8180

# Launch Shiny app
CMD ["R", "-e", "library('shiny') ; options(shiny.port = 8180, shiny.host = '0.0.0.0') ; shiny::runApp('Shiny_Acomys')"]
