#Set rotational / angular momenta from selected rotational state
#NEEDS WORK
function calc_ROTN(AM,EROT)
   #Prior name = ROTN

   # C         CALCULATE ANGULAR MOMENTUM, MOMENT OF INERTIA TENSOR,
   # C         ANGULAR VELOCITY, AND ROTATIONAL ENERGY

   # C         CALCULATE ANGULAR MOMENTUM.  THE CENTER OF MASS COORDINATES
   # C         QQ AND MOMENTA PP COME FROM SUBROUTINE CENMAS THROUGH COMMON
   # C         BLOCK WASTE.
   if (NAM.EQ.0) THEN
      AM(1)=0.0D0
      AM(2)=0.0D0
      AM(3)=0.0D0
      AM(4)=0.0D0
      EROT=0.0D0

      do I=1,2
         J=L(I)
         J3=3*J
         J2=J3-1
         J1=J2-1
         AM(1)=AM(1)+(QQ(J2)*PP(J3)-QQ(J3)*PP(J2))
         AM(2)=AM(2)+(QQ(J3)*PP(J1)-QQ(J1)*PP(J3))
         AM(3)=AM(3)+(QQ(J1)*PP(J2)-QQ(J2)*PP(J1))
      end
      AM(4)=SQRT(AM(1)**2+AM(2)**2+AM(3)**2)
      
      AIXX=0.0D0
      do I=1,2
         J=3*L(I)+1
         SR=0.0D0
         do K=1,3
            SR=SR+QQ(J-K)**2
         end
         AIXX=AIXX+SR*W(L(I))
      end
      EROT=AM(4)**2/AIXX/2.0D0/C1
      return
   end
   #         CALCULATE THE MOMENT OF INERTIA TENSOR
   AIXX=0.0D0
   AIYY=0.0D0
   AIZZ=0.0D0
   AIXY=0.0D0
   AIXZ=0.0D0
   AIYZ=0.0D0
   do I=1,2
      J=L(I)
      J3=3*J
      J2=J3-1
      J1=J2-1
      AIXX=AIXX+W(J)*(QQ(J2)**2+QQ(J3)**2)
      AIYY=AIYY+W(J)*(QQ(J1)**2+QQ(J3)**2)
      AIZZ=AIZZ+W(J)*(QQ(J1)**2+QQ(J2)**2)
      AIXY=AIXY+W(J)*QQ(J1)*QQ(J2)
      AIXZ=AIXZ+W(J)*QQ(J1)*QQ(J3)
      AIYZ=AIYZ+W(J)*QQ(J2)*QQ(J3)
   end
   DET=AIXX*(AIYY*AIZZ-AIYZ*AIYZ)-AIXY*(AIXY*AIZZ+AIYZ*AIXZ)-
   *AIXZ*(AIXY*AIYZ+AIYY*AIXZ)
   #
   #         CALCULATE INVERSE OF THE INERTIA TENSOR
   #
   if (ABS(DET).GE.0.01D0) THEN
      UXX=(AIYY*AIZZ-AIYZ*AIYZ)/DET
      UXY=(AIXY*AIZZ+AIXZ*AIYZ)/DET
      UXZ=(AIXY*AIYZ+AIXZ*AIYY)/DET
      UYY=(AIXX*AIZZ-AIXZ*AIXZ)/DET
      UYZ=(AIXX*AIYZ+AIXZ*AIXY)/DET
      UZZ=(AIXX*AIYY-AIXY*AIXY)/DET
      #     CALCULATE ANGULAR VELOCITIES
      WX=UXX*AM(1)+UXY*AM(2)+UXZ*AM(3)
      WY=UXY*AM(1)+UYY*AM(2)+UYZ*AM(3)
      WZ=UXZ*AM(1)+UYZ*AM(2)+UZZ*AM(3)
   else
      #    CALCULATE ROTATIONAL ENERGY
      AIXX=0.0D0
      do I=1,2
         J=3*L(I)+1
         SR=0.0D0
         do K=1,3
               SR=SR+QQ(J-K)**2
         end
         AIXX=AIXX+SR*W(L(I))
      end
      EROT=AM(4)**2/AIXX/2.0D0/C1
      RETURN
   end
   EROT=(WX*AM(1)+WY*AM(2)+WZ*AM(3))/2.0D0/C1
   RETURN
end