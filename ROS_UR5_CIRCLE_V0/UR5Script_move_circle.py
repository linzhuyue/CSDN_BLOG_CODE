#!/usr/bin/env python
# -*- coding: utf-8 -*-
import rospy
import numpy
from std_msgs.msg import String
from frompitoangle import *
from ur5_kinematics import *
from ur5_pose_get import *
from transfer import *

class UrCircle:
    def __init__(self,weights,radius):
        self.weights = weights
        self.radius=radius#m
        self.cont=150
    def Init_node(self):
        rospy.init_node("move_ur5_circle")
        pub = rospy.Publisher("/ur_driver/URScript", String, queue_size=10)
        return pub
    def get_urobject_ur5kinetmatics(self):
        ur0 = Kinematic()
        return ur0
    def get_draw_circle_xy(self,t,xy_center_pos):
        x = xy_center_pos[0] + self.radius * math.cos( 2 * math.pi * t / self.cont )
        y = xy_center_pos[1] + self.radius * math.sin( 2 * math.pi * t / self.cont)
        return  [x,y]
    def get_IK_from_T(self,T,q_last):
        ur0 = self.get_urobject_ur5kinetmatics()
        return ur0.best_sol(self.weights,q_last,T)
    def get_q_list(self,T_list,qzero):
        ur0 = self.get_urobject_ur5kinetmatics()
        tempq=[]
        resultq=[]
        for i in xrange(self.cont):
            if i==0:
                tempq=qzero
                firstq = ur0.best_sol(self.weights, tempq, T_list[i])
                tempq=firstq
                resultq.append(firstq.tolist())
                print "firstq", firstq
            else:
                qq = ur0.best_sol(self.weights, tempq, T_list[i])
                tempq=qq
                print "num i qq",i,qq
                resultq.append(tempq.tolist())
        return resultq
    # T is numpy.array
    def get_T_translation(self, T):
        trans_x = T[3]
        trans_y = T[7]
        trans_z = T[11]
        return [trans_x, trans_y, trans_z]
    def insert_new_xy(self,T,nx,ny):
        temp=[]
        for i in xrange(12):
            if i==3:
                temp.append(nx)
            elif i==7:
                temp.append(ny)
            else:
                temp.append(T[i])
        return temp
    def get_new_T(self,InitiT,xy_center_pos):
        temp=[]
        for i in xrange(self.cont):
            new_xy=self.get_draw_circle_xy(i,xy_center_pos)
            new_T=self.insert_new_xy(InitiT,new_xy[0],new_xy[1])
            #print "initial T\n",InitiT
            #print "new_T\n",i,new_T
            temp.append(new_T)
        return temp
    def urscript_pub(self, pub, qq, vel, ace, t):

        ss = "movej([" + str(qq[0]) + "," + str(qq[1]) + "," + str(qq[2]) + "," + str(
            qq[3]) + "," + str(qq[4]) + "," + str(qq[5]) + "]," + "a=" + str(ace) + "," + "v=" + str(
            vel) + "," + "t=" + str(t) + ")"
        print("---------ss:", ss)
            # ss="movej([-0.09577000000000001, -1.7111255555555556, 0.7485411111111111, 0.9948566666666667, 1.330836666666667, 2.3684322222222223], a=1.0, v=1.0,t=5)"
        pub.publish(ss)
def main():
    t=0
    vel=0.2
    ace=50
    qstart=[133.70,-79.09,103.45,-111.83,-89.54,-212.29]
    # qq=[45.91,-72.37,61.52,-78.56,-90.49,33.71]
    # q=display(getpi(qq))
    ratet = 30
    radius=0.3
    weights = [1.] * 6
    T_list=[]
    urc=UrCircle(weights,radius)
    pub=urc.Init_node()
    rate = rospy.Rate(ratet)

    # first step go to initial pos
    q = display(getpi(qstart))
    # urc.urscript_pub(pub,q,vel,ace,t)
    # second get T use urkinematics
    urk = urc.get_urobject_ur5kinetmatics()
    F_T = urk.Forward(q)
    print "F_T", F_T
    TransT = urc.get_T_translation(F_T)
    print "TransT", TransT
    xy_center_pos = [TransT[0], TransT[1]]
    print "xy_center_pos", xy_center_pos
    T_list = urc.get_new_T(F_T, xy_center_pos)
    print "T_list", T_list
    reslut_q = urc.get_q_list(T_list, q)
    print "reslut_q", reslut_q
    cn=0
    while not rospy.is_shutdown():
        urc.urscript_pub(pub, reslut_q[cn], vel, ace, t)
        cn+=1
        if cn==len(reslut_q):
            cn=0
        # print "cn-----\n",reslut_q[cn]

        rate.sleep()
if __name__ == '__main__':
        main()