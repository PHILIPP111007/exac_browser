# Usage


*If you would like to use the ExAC browser, the most recent stable version is hosted at <http://exac.broadinstitute.org>*

Advanced: The following instructions are useful for cloning the browser (e.g. to load custom sites/coverage data).
Most users will not need to go through this process.

# Installation

First (as this can run in parallel), get the datasets that the browser uses and put them into an 'exac_data' directory:

    wget http://broadinstitute.org/~konradk/exac_browser/exac_browser.tar.gz .
    tar zxvf exac_browser.tar.gz
    cd ..

## Dependencies

Follow these instructions to get Python and Homebrew installed on your Mac:
<http://docs.python-guide.org/en/latest/starting/install/osx/>

Install MongoDB:

```sh
brew install mongodb
# or
sudo port install mongodb
```

Create a directory to hold your mongo database files:

```sh
mkdir database
```

In a separate tab, start the mongo database server:

```sh
mongod --dbpath database
# mongod --port 27017 --dbpath ./database --replSet myreplicaset
```

This local server needs to be running at all times when you are working on the site.
You could do this in the background if you want or set up some startup service,
but I think it's easier just to open a tab you can monitor.

```sh
micromamba env create -f ./env.yml
```

## Setup

At this point, it's probably worth quickly checking out the code structure if you haven't already :)

Now we must load the database from those flat files.
This is a single command, but it can take a while (can take advantage of parallel loads by modifying LOAD\_DB\_PARALLEL\_PROCESSES in exac.py):

    python manage.py load_db

You won't have to run this often - most changes won't require rebuilding the database.
That said, this is (and will remain) idempotent,
so you can run it again at any time if you think something might be wrong - it will reload the database from scratch.
You can also reload parts of the database using any of the following commands:

```sh
python manage.py load_variants_file
python manage.py load_dbsnp_file
python manage.py load_base_coverage
python manage.py load_gene_models
```

Then run:

```sh
python manage.py precalculate_metrics
```

Then, you need to create a cache for autocomplete and large gene purposes:

```sh
python manage.py create_cache
```

## Running the site

Note that if you are revisiting the site after a break, make sure your virtualenv is `activate`'d.

You can run the development server with:

```sh
python main.py
```

And visit on your browser:

    http://localhost:8000
    http://localhost:8000/gene/ENSG00000237683
    http://localhost:8000/variant/20-76735-A-T

For testing, you can open up an interactive shell with:

```sh
python manage.py shell
```
