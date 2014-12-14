Carbonate 2.0
=============

Installation
------------

Clone the repository:

    git clone https://www.github.com/csvoss/carbonate
    git config --global user.name "Your Name"
    git config --global user.email youremail@example.com

Open a command line in the same folder as the repository, and start Vagrant. (Problems? Try reinstalling VirtualBox.) (No command Vagrant found? If you're using a Mac, try running: `export PATH=/Applications/Vagrant/bin:$PATH`.)

    vagrant up
    vagrant ssh

Now you're in the virtual machine! It has Django installed, but you need to install some other things. Open requirements.txt and run the commands there. For example:

    sudo apt-get install python-openbabel
    sudo apt-get install sqlite3
    sudo apt-get install python-django-south
    sudo pip install django-forkit
    sudo pip install django-crispy-forms
    sudo pip install rply

To test that OpenBabel installed correctly:

    python
    >>> import openbabel

Now, try to run the Django server. This code will run the server locally, at port 4000.

    cd /vagrant
    cd orgo_django/
    python manage.py runserver 0.0.0.0:4000

If everything goes well, you will be able to access your copy running on the local server. Open `localhost:4000` in your browser.

Project Layout
==============

This Django project is divided into three folders: /api, /app, and /orgo.

API
---
This contains the chemistry engine of Carbonate, which other parts of the website can access through its API. It knows about reagents, reactions, and molecules. It can convert molecules to and from the SMILES representation, and export them as SVG images.

App
---
This contains standalone user interfaces -- for example, the Reaction Tutor, simple single-step problems, and predict-the-products problems -- that could easily be embedded as iframes in other pages.

Orgo
----
This contains the code of the orgo.mit.edu site, mostly unmodified except for those parts that were refactored into the API.
