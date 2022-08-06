# Install the X server

This process is necessary to get the GUI (Graphical User Interface) to work properly

Follow the following step:

1. Go to the tools drop down menu and select install software.
2. Scroll down and select "X server" and then click the little download icon all the way to the right. 
3. Let that finish (which could take up to 10 minutes) and then restart the box. (Go to the Project drop down and select "restart box")

At this point when you start the virtual desktop as indicated in the main README, you should see a black screen.  When you run the Display class from the terminal, your GUI will show up on this Screen.

## Supplemental Reset process

If you terminate the installation process early subsequent installation attempts will fail. Follow this processes to reset things:

1. In the Project drop down menu, select "Reset Box." Then type in the number they display to confirm. This will reset the underlying virtual machine.
2. Once that is complete, go back and do the X server installation per the previous instructions. Make sure to leave the window open and wait until it says you can restart your box. This could take as long as 10 minutes.
3. Once it fully finishes, go back to the Project down menu and select "Restart Box." Once that's done you should be good to go.