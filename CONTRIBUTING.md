# Contributing to *ASLPrep*

Welcome to the *ASLPrep* repository!
We're excited you're here and want to contribute.

**Imposter's syndrome disclaimer**[^1]: We want your help. No, really.

There may be a little voice inside your head that is telling you that
you're not ready to be an open-source contributor; that your skills
aren't nearly good enough to contribute. What could you possibly offer a
project like this one?

We assure you - the little voice in your head is wrong. If you can
write code at all, you can contribute code to open-source. Contributing
to open-source projects is a fantastic way to advance one's coding
skills. Writing perfect code isn't the measure of a good developer (that
would disqualify all of us!); it's trying to create something, making
mistakes, and learning from those mistakes. That's how we all improve,
and we are happy to help others learn.

Being an open-source contributor doesn't just mean writing code, either.
You can help out by writing documentation, tests, or even giving
feedback about the project (and yes - that includes giving feedback
about the contribution process). Some of these contributions may be the
most valuable to the project as a whole, because you're coming to the
project with fresh eyes, so you can see the errors and assumptions that
seasoned contributors have glossed over.

  1. The tool only and fully supports BIDS and BIDS-Derivatives for the input and output data.
  1. The tool is packaged as a fully-compliant [BIDS-App](https://bids-apps.neuroimaging.io), not just in its user
     interface, but also in the continuous integration, testing and delivery.
  1. The ASL is not fully supported by BIDS as of now, it is in [progress](https://docs.google.com/document/d/15tnn5F10KpgHypaQJNNGiNKsni9035GtDqJzWqkkP6c/edit#). A simple tool has been developed to convert ASL to BIDS format, [aslbids](https://github.com/PennLINC/aslbids).
  1. The tool is rigorously restricted to ASL preprocessing and cbf computation, including (but not limited to):
     input/output metadata assessment, head-motion correction, susceptibility-distortion correction,
     co-registration with anatomical data, spatial normalization to neuroimaging templates.
     In other words, the tool does not deal with time filtering*, smoothing*, modeling,
     or connectivity extraction.
  1. The tool is **agnostic to subsequent analysis**, i.e., any software supporting BIDS-Derivatives
     for its inputs should be amenable to fit GLMs, extract time-series for connectivity analyses, etc.
  1. The tool is thoroughly and transparently documented (including the generation of individual reports
     that can be used as scaffolds for understanding the underpinnings and design decisions of the tool).
  1. The tool is community-driven, with a very open concept of contribution that is always credited
     with authorship offers when writing relevant papers.

## Practical guide to submitting your contribution

These guidelines are designed to make it as easy as possible to get involved.
If you have any questions that aren't discussed below,
please let us know by opening an [issue](https://github.com/PennLINC/aslprep/issues)!

Before you start, you'll need to set up a free [GitHub](https://github.com) account and sign in.
Here are some [instructions](https://docs.github.com/en/github/getting-started-with-github/signing-up-for-a-new-github-account).

Already know what you're looking for in this guide? Jump to the following sections:

* [Joining the conversation](#joining-the-conversation)
* [Contributing through Github](#contributing-through-github)
* [Understanding issues](#understanding-issues)
* [Making a change](#making-a-change)
* [Structuring contributions](#ASLPrep-coding-style-guide)
* [Licensing](#licensing)
* [Recognizing contributors](#recognizing-contributions)

## Joining the conversation

*ASLPrep* is maintained by a growing group of enthusiastic developers&mdash;
and we're excited to have you join!
Most of our discussions will take place on open [issues](https://github.com/PennLINC/aslprep/issues).

We also encourage users to report any difficulties they encounter on [NeuroStars](https://neurostars.org),
a community platform for discussing neuroimaging.

We actively monitor both spaces and look forward to hearing from you in either venue!

## Contributing through GitHub

[git][
_git] is a really useful tool for version control.
[GitHub][link_github] sits on top of git and supports collaborative and distributed working.

If you're not yet familiar with `git`, there are lots of great resources to help you *git* started!
Some of our favorites include the [git Handbook](https://guides.github.com/introduction/git-handbook/) and
the [Software Carpentry introduction to git](https://swcarpentry.github.io/git-novice/).

On GitHub, You'll use [Markdown](https://www.markdownguide.org) to chat in issues and pull requests.
You can think of Markdown as a few little symbols around your text that will allow GitHub
to render the text with a little bit of formatting.
For example, you could write words as bold (`**bold**`), or in italics (`*italics*`),
or as a [link](https://youtu.be/dQw4w9WgXcQ) (`[link](https://youtu.be/dQw4w9WgXcQ)`) to another webpage.

GitHub has a really helpful page for getting started with
[writing and formatting Markdown on GitHub](https://docs.github.com/en/github/writing-on-github/basic-writing-and-formatting-syntax).

## Understanding issues

Every project on GitHub uses [issues](https://github.com/PennLINC/aslprep/issues) slightly differently.

The following outlines how the *ASLPrep* developers think about these tools.

* **Issues** are individual pieces of work that need to be completed to move the project forward.
A general guideline: if you find yourself tempted to write a great big issue that
is difficult to describe as one unit of work, please consider splitting it into two or more issues.

    Issues are assigned [labels](#issue-labels) which explain how they relate to the overall project's
    goals and immediate next steps.

### Issue Labels

The current list of issue labels are [here](https://github.com/PennLINC/aslprep/labels) and include:

* These issues contain a task that is amenable to new contributors because it doesn't entail a steep learning curve.*

    If you feel that you can contribute to one of these issues,
    we especially encourage you to do so!

    If you find new a bug, please give as much detail as possible in your issue,
    including steps to recreate the error.
    If you experience the same bug as one already listed,
    please add any additional information that you have as a comment.

 *These issues are asking for new features and improvements to be considered by the project.*

    Please try to make sure that your requested feature is distinct from any others
    that have already been requested or implemented.
    If you find one that's similar but there are subtle differences,
    please reference the other request in your issue.

In order to define priorities and directions in the development roadmap,
we have two sets of special labels:


## Making a change

We appreciate all contributions to *ASLPrep*,
but those accepted fastest will follow a workflow similar to the following:

1. **Comment on an existing issue or open a new issue referencing your addition.**<br />
  This allows other members of the *ASLPrep* development team to confirm that you aren't
  overlapping with work that's currently underway and that everyone is on the same page
  with the goal of the work you're going to carry out.<br />
  [This blog](http://alblue.bandlem.com/2011/03/git-tip-of-week-pushing-and-pulling.html) is a nice explanation of why putting this work in up front
  is so useful to everyone involved.

## Recognizing contributions

We welcome and recognize all contributions regardless their size, content or scope:
from documentation to testing and code development.
You can see a list of current developers and contributors in our [zenodo file](https://zenodo.org/record/3932697#.XxeLpC2z1wc).
Before every release, a new [zenodo file](https://zenodo.org/record/3932697#.XxeLpC2z1wc) will be generated.
The [update script][link_update_script] will also sort creators and contributors by
the relative size of their contributions, as provided by the `git-line-summary` utility
distributed with the `git-extras` package.
Last positions in both the *creators* and *contributors* list will be reserved to
the project leaders.
These special positions can be revised to add names by punctual request and revised for
removal and update of ordering in an scheduled manner every two years.
All the authors enlisted as *creators* participate in the revision of modifications.

### Developers

Developers are members of a wonderful team _driving the project_.
Names and contacts of all developers are included in the
Examples of steering activities that _drive the project_ are: actively participating in the
follow-up meetings, leading documentation sprints, helping in the design of the tool and definition of the roadmap,
providing resources (in the broad sense, including funding), code-review, etc.

### Contributors

Contributors are people that are 
actively help or have previously helped the project in a broad sense: writing code, writing documentation,
benchmarking modules of the tool, proposing new features, helping improve the scientific
rigor of implementations, giving out support on the different communication
channels and GitHub issues etc..
If you are new to the project, don't forget to add your name and affiliation to the list
of contributors there!
Before every release, unlisted contributors will be invited again to add their names to the file
(just in case they missed the automated message from our Welcome Bot).

### Publications

Anyone listed as a *developer* or a *contributor* can start the submission process of a manuscript
as first author.
To compose the author list, all the *creators* MUST be included (except for those people who
opt to drop-out) and all the *contributors* MUST be invited to participate.
First authorship(s) is (are) reserved for the authors that originated and kept the initiative
of submission and wrote the manuscript.
To generate the ordering of your paper, please run ``python .maint/paper_author_list.py`` from the
root of the repository, on the up-to-date ``upstream/master`` branch.
Then, please modify this list and place your name first.
All developers and contributors are pulled together in a unique list, and last authorships assigned.
*fMRIPrep* and its community adheres to open science principles, such that a pre-print should
be posted on an adequate archive service (e.g., [ArXiv](https://arxiv.org) or
[BioRxiv](https://biorxiv.org)) prior publication.


## Licensing

*ASLPrep* is licensed under the Apache 2.0 license.
By contributing to *ASLPrep*,
you acknowledge that any contributions will be licensed under the same terms.

## Thank you!

You're awesome. :wave::smiley:
[Github](https://github.com/): https://github.com/
[ASLPrep](https://github.com/pennlinc/aslprep): https://github.com/pennlinc/aslprep
