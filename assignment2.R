# check your set-up

usethis::git_remotes()
# you should see an origin (in your github account) and upstream (the other's github account)

usethis::git_sitrep()

# _ Code your own solution -----

# you will do this everytime you have a chunck of code that you want to sink with your GitHub repository online. I keep this GitHub_software_tools.R file in a central location, and always open it after I have openend the project file.

# CREATE A NEW BRANCH - do this after every pull request has been merged!

usethis::pr_init(branch = "lvillasa") # change formidable w/ more informative branch name!

usethis::git_sitrep()
# check that Git local (project) Current local branch is lvillasa
# and this is differet from the git_sitrep output from line38!

# now do your R work
print("second test")
# for this assignment, create a script file named
# <your-university-user-name-from-your-uog-email-address>.R
# for me, that would be cottenie.R
# this file will consist of your commented script file for the randomforest code
# put it in the correct folder!

# COMMIT your new file with the changes in the "git" tab next to Environmnet

# _ Submit pull request -----

usethis::pr_push()

# in the launched browser window, click "Create pull request"
# in the description, add - "This is my random forest script file :-)"
# click "Create pull request"

# then you wait for me to review and accept - this is your "quiz" for this week

# after I have merged the pull request

usethis::pr_pull() # make sure your branch is up-to-date

#usethis::pr_finish() # cleaning everything up, deleting local 

#usethis::git_sitrep()
