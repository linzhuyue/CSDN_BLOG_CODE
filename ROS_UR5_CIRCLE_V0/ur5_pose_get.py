#!/usr/bin/env python
import rospy
from std_msgs.msg import String
# from ar_track_alvar_msgs.msg import AlvarMarkers
#from std_msgs.msg import Empty
from sensor_msgs.msg import JointState
from frompitoangle import *
import os, time
import sys

class Urposition():

    def __init__(self, name = "ur_info_subscriber" ):

        self.name = name
        # self.pos_dict = {}
        self.ur_pose_buff_list = []
        self.ave_ur_pose = []
        self.now_ur_pos = []
        self.tmp_sum = [0]*6

        # pass

    def Init_node(self):
        rospy.init_node(self.name)
        sub = rospy.Subscriber("/joint_states", JointState, self.callback)
        return sub

    def urscript_pub(self, pub, qq, vel, ace, t, rate):

        rate = rospy.Rate(rate)

        ss = "movej([" + str(qq[0]) + "," + str(qq[1]) + "," + str(qq[2]) + "," + str(
            qq[3]) + "," + str(qq[4]) + "," + str(qq[5]) + "]," + "a=" + str(ace) + "," + "v=" + str(
            vel) + "," + "t=" + str(t) + ")"
        print("---------ss:", ss)
            # ss="movej([-0.09577000000000001, -1.7111255555555556, 0.7485411111111111, 0.9948566666666667, 1.330836666666667, 2.3684322222222223], a=1.0, v=1.0,t=5)"
        pub.publish(ss)

        rate.sleep()
                # rostopic pub /ur_driver/URScript std_msgs/String "movej([-0.09577000000000001, -1.7111255555555556, 0.7485411111111111, 0.9948566666666667, 1.330836666666667, 2.3684322222222223], a=1.4, v=1.0)"
        # time.sleep(3)



    def callback(self, msg):
        self.read_pos_from_ur_joint( msg )
        self.ur_pose_buff_list, self.ave_ur_pose = self.pos_filter_ur( self.ur_pose_buff_list, self.now_ur_pos )


    def read_pos_from_ur_joint(self, msg):
        self.now_ur_pos = list( msg.position )
        # print("pos_msg:", pos_msg)
        # ur_pose = [ pos_msg.x , pos_msg.y, pos_msg.z, quaternion_msg.x, quaternion_msg.y, quaternion_msg.w  ]
        # return pos_msg

    def pos_filter_ur(self, pos_buff, new_data ):
        # tmp_sum = [0]*6
        ave_ur_pose = [0]*6
        if len( pos_buff ) == 10 :
            # print("new_data:", new_data)
            # print("tmp_sum before:", tmp_sum)
            # print("pos_buff[0]:", pos_buff[0])
            res1 = list_element_minus( self.tmp_sum , pos_buff[0] )
            self.tmp_sum = list_element_plus( res1 , new_data )
            # print("----------res1:", res1)
            # print("----------tmp_sum after:", tmp_sum)
            pos_buff = pos_buff[1:]
            pos_buff.append(new_data)
            ave_ur_pose = list_element_multiple( self.tmp_sum, 1.0/10 )
            # print( "len:", len( pos_buff ))
            # print ("10----ave_pos_ur:", ave_ur_pose)
            # time.sleep(2)
            # sys.exit(0)
        else:
            pos_buff.append(new_data)
            self.tmp_sum = list_element_plus( self.tmp_sum, new_data )
            ave_ur_pose = pos_buff[-1] # get the last element
            # print ("----------tmp_sum:", self.tmp_sum)
        # pos_buff.append( new_data )
        # print("---------------len:", pos_buff)
        return pos_buff, ave_ur_pose



def list_element_plus( v1, v2):
    res = list(map( lambda x: x[0] + x[1] , zip(v1,v2)))
    # print "plus :", res
    return res

def list_element_minus( v1, v2):
    res = list(map( lambda x: x[0] - x[1] , zip(v1,v2)))
    # print "minus :", res
    return res

def list_element_multiple( v1, num ):
    return [ item * num  for item in v1 ]


def main():
    ur_info_reader = Urposition()
    ur_info_reader.Init_node()
    while not rospy.is_shutdown():
        pass
        # print ( "now_pos: ", type(ur_info_reader.now_ur_pos))

        # print ("ave_pos_ur:", ur_info_reader.ave_ur_pose)
    # rospy.spin()


if __name__ == "__main__":
    # v1 = [1,2,43.0,5]
    # v2 = [3, 1, 0.9, 2.9 ]
    # print list_element_minus( v1 , v2 )
    main()