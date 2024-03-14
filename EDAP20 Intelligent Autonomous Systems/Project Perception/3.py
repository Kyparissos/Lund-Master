#! /usr/bin/env python

# Import the Python library for ROS
import rospy
# Import the Int32 message from the std_msgs package
from std_msgs.msg import Int32    

class Twist:
    def __init__(move,linear,angular):
        move.linear.x = linear
        move.angular.z = angular
cmd = Twist(0.5,1.0)

# Initiate a Node named 'topic_publisher'
rospy.init_node('topic_publisher')

# Create a Publisher object, that will publish on the /counter topic
# messages of type Int32

pub = rospy.Publisher('/cmd_vel ', Twist , queue_size=1)    
                                           
# Set a publish rate of 2 Hz
#rate = rospy.Rate(2)
# Create a variable of type Int32
#count = Int32()
# Initialize 'count' variable
#count.data = 0                    

# Create a loop that will go until someone stops the program execution
while not rospy.is_shutdown():
  # Publish the message within the 'count' variable
  pub.publish(cmd)
  # Increment 'count' variable
  #count.data += 1
  # Make sure the publish rate maintains at 2 Hz
  rate.sleep()                             