====
Git
====

---------------------
Setting up git config
---------------------

.. code-block:: console

    git config --global user.name
    git config --global user.email
    ...

--------------
Basic Workflow
--------------


.. image:: ../images/git_basics.png

.. code-block:: console

    git init
    $ git add
    $ git commit -m
    $ git push
    $ git pull

---------
Branching
---------

.. code-block:: console

    git status
    git branch
    git checkout -b issue_1
    git fetch
    git rebase master
    git checkout master
    git merge
    git branch -d issue_1

--------------------------
Working through an example
--------------------------
