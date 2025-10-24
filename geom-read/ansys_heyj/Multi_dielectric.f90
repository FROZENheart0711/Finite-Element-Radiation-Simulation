         subroutine readcdb(tmp_char)
         implicit none
	 integer :: maxnode,maxpatch,maxfacepatch,maxlv
         character*256 :: Path,Dataname
		 integer :: Uinp=1,Uout_xyz=2,Uout_ipat=3,Uout_lv=4

         character*256 :: tmp_filename1,tmp_filename2
		 integer :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
		 
		 integer :: i,j,k,ii,jj,kk
		 real :: temp1,temp2,temp3,temp4
		 complex :: cmp_tmp1,cmp_tmp2,cmp_tmp3,cmp_tmp4
		 character*256 :: tmp_char1,tmp_char2,tmp_char3
                 character(len=*)::tmp_char
         integer nargs, iargc

	 integer,allocatable,dimension(:) :: arr1_tmp1,arr1_tmp2,arr1_tmp3,&
		                                     arr1_tmp4,arr1_tmp
		 integer,allocatable,dimension(:,:) :: arr2_tmp1,arr2_tmp2

		 real,allocatable,dimension(:) :: arr_temp1,arr_temp2,arr_temp3
		 integer,allocatable,dimension(:) :: tp1,tp2,tp3,tp4

		 
		 !c.....get commend line 

		 write(*,*) "transforming the *.cdb to *.xyz,*.ipat,*.lv";
		 write(*,*) "waiting....."
		 tmp_filename2='.cdb'
                 call fstrcat (tmp_char,tmp_filename2,tmp_filename1)
		 open( Uinp,file=tmp_filename1 )
                 print*,'file read ok'
		       do ii=1,5
                  read(Uinp,*) 
			   enddo
			   read(Uinp,*) tmp_char1,tmp_char2,tmp1
			   
			   if (tmp1>0 ) then
			       maxnode=tmp1
			   else
			       write(*,*) "Maxnode is wrong at line 6."
				   stop
			   endif

			   read(Uinp,*) tmp_char1,tmp_char2,tmp1
			   
			   if (tmp1>0 ) then
			       maxpatch=tmp1
			   else
			       write(*,*) "Maxpatch is wrong at line 7."
				   stop
			   endif
               
			   do ii=1,2
                  read(Uinp,*) 
			   enddo

			   read(Uinp,*) tmp_char1,tmp_char2

			   tmp_char3=trim(tmp_char1)//','//trim(tmp_char2);
			   if (trim(tmp_char3)/='(3i8,6e16.9)') then
			       write(*,*) " The file's format dosen't conform..."
				   write(*,*) " The program's format is (3i8,6e16.9), but &
				                the format of the *.cdb (line 10) is: ",trim(tmp_char3)
		write(*,*) " Please check the *.cdb again or change the format &
			                of this transform program "
		       stop
			   endif
               
			   

               
			   !============================================
			   !  Read and transform the *.xyz
			   !============================================
			   tmp1=1;   tmp2=maxnode;
			   allocate( arr_temp1(tmp1:tmp2),arr_temp2(tmp1:tmp2),&
						 arr_temp3(tmp1:tmp2) );
			   
			   arr_temp1=0; arr_temp2=0; arr_temp3=0;

			   do ii=1,maxnode
			      read(Uinp,1201) tmp1,tmp2,tmp3,&
				                  arr_temp1(ii),arr_temp2(ii),arr_temp3(ii) 
               enddo
                                   tmp_filename2='.xyz'
                 call fstrcat (tmp_char,tmp_filename2,tmp_filename1)
			   open(Uout_xyz,file=tmp_filename1)
			        write(Uout_xyz,*) maxnode
					do ii=1,maxnode
					   write(Uout_xyz,*) arr_temp1(ii),arr_temp2(ii),arr_temp3(ii)
					enddo
			   close(Uout_xyz);
               deallocate( arr_temp1,arr_temp2,arr_temp3 );
               
			   write(*,*) "finish creating ",trim(tmp_filename1);
			 
			   !============================================
			   !  Read the *.ipat
			   !============================================
               
			   
			   read(Uinp,*)
			   read(Uinp,*)
			   read(Uinp,*) tmp_char1

			   if (trim(tmp_char1)/='(19i8)') then
			       write(*,*) " The file's format dosen't conform..."
				   write(*,*) " The program's format is (19i8), but &
				                the format of the *.cdb (line 849) is: ",trim(tmp_char1)
				   write(*,*) " Please check the *.cdb again or change the format &
				                of this transform program "
			       stop
			   endif
               
			   tmp1=1;   tmp2=maxpatch;
			   allocate( arr1_tmp1(tmp1:tmp2),arr1_tmp2(tmp1:tmp2),&
						 arr1_tmp3(tmp1:tmp2) )
			   arr1_tmp1=0; arr1_tmp2=0; arr1_tmp3=0;

			   allocate( arr1_tmp(1:19) ); arr1_tmp=0;
               
			   do ii=1,maxpatch
			      tmp1=arr1_tmp(2);
			      read(Uinp,1203) arr1_tmp;
				  tmp2=arr1_tmp(2);
				  arr1_tmp1(ii)=arr1_tmp(12);
				  arr1_tmp2(ii)=arr1_tmp(13);
				  arr1_tmp3(ii)=arr1_tmp(14); 
				  if ( ii>=3 .and. tmp1/=tmp2 ) then
				       maxfacepatch=ii-1;
					   if ( ii/=arr1_tmp(11) ) then
					        write(*,*) "Sth. wrong about the *.cdb,&
							            the maximum facepatch is wrong,&
										please check the *.cdb file or check"
					        stop
					   endif

				       exit
			      endif
               enddo
			   
               
			   allocate( arr2_tmp1(1:maxfacepatch,1:4) ); arr2_tmp1=0;
			   arr2_tmp1(:,1)=arr1_tmp1;
			   arr2_tmp1(:,2)=arr1_tmp2;
			   arr2_tmp1(:,3)=arr1_tmp3;
               deallocate( arr1_tmp1,arr1_tmp2,arr1_tmp3);
               
               
			   !============================================
			   !  Read the *.lv
			   !============================================
               maxlv=maxpatch-maxfacepatch;

			   allocate( arr2_tmp2(1:maxlv,1:9) ); arr2_tmp2=0;

			   do ii=1,maxlv

			      arr2_tmp2(ii,1:4)=arr1_tmp(12:15);
				  arr2_tmp2(ii,5)=arr1_tmp(2);
			      read(Uinp,1203) arr1_tmp;
               enddo
               close(Uinp);  deallocate( arr1_tmp );
               

			   !============================================
			   !  Judge the common triangle face.....
			   !  By four operation
			   !
			   !  x+y+z;
			   !  x*y*z;
			   !  x**2+y**2+z**2
			   !  x**3+y**3+z**3 
			   !============================================

			   allocate( arr1_tmp1(1:maxfacepatch),arr1_tmp2(1:maxfacepatch),&
			             arr1_tmp3(1:maxfacepatch),arr1_tmp4(1:maxfacepatch) );
			   
			   arr1_tmp1=0; arr1_tmp2=0; arr1_tmp3=0; arr1_tmp4=0;

			   do ii=1,maxfacepatch

                  tmp1=arr2_tmp1(ii,1);
				  tmp2=arr2_tmp1(ii,2);
				  tmp3=arr2_tmp1(ii,3);
                  
				  Call Operate_four(tmp1,tmp2,tmp3,&
                                    arr1_tmp1(ii),arr1_tmp2(ii),&
									arr1_tmp3(ii),arr1_tmp4(ii))
		
               enddo

			   allocate( tp1(1:4),tp2(1:4),tp3(1:4),tp4(1:4) ); 
			   tp1=0; tp2=0; tp3=0; tp4=0;
               
			   do jj=1,maxlv

                  tmp1=arr2_tmp2(jj,2)
				  tmp2=arr2_tmp2(jj,3)
				  tmp3=arr2_tmp2(jj,4)
                  
				  Call Operate_four(tmp1,tmp2,tmp3,&
                                    tp1(1),tp2(1),tp3(1),tp4(1))
                  
                  tmp1=arr2_tmp2(jj,1)
				  tmp2=arr2_tmp2(jj,3)
				  tmp3=arr2_tmp2(jj,4)
                  
				  Call Operate_four(tmp1,tmp2,tmp3,&
                                    tp1(2),tp2(2),tp3(2),tp4(2))

				  tmp1=arr2_tmp2(jj,2)
				  tmp2=arr2_tmp2(jj,1)
				  tmp3=arr2_tmp2(jj,4)
                  
				  Call Operate_four(tmp1,tmp2,tmp3,&
                                    tp1(3),tp2(3),tp3(3),tp4(3))
                  
				  tmp1=arr2_tmp2(jj,2)
				  tmp2=arr2_tmp2(jj,3)
				  tmp3=arr2_tmp2(jj,1)
                  
				  Call Operate_four(tmp1,tmp2,tmp3,&
                                    tp1(4),tp2(4),tp3(4),tp4(4))


			   do ii=1,maxfacepatch
			      
				  tmp1=arr1_tmp1(ii);
				  tmp2=arr1_tmp2(ii);
				  tmp3=arr1_tmp3(ii);
				  tmp4=arr1_tmp4(ii);

                  if (tmp1==tp1(1)) then
                      if (tmp2==tp2(1)) then
                          if (tmp3==tp3(1)) then
                              if (tmp4==tp4(1)) then
                                  arr2_tmp2(jj,6)=1;
								  arr2_tmp1(ii,4)=1;
							  endif
                          endif
                      endif
                  endif
                  if (tmp1==tp1(2)) then
                      if (tmp2==tp2(2)) then
                          if (tmp3==tp3(2)) then
                              if (tmp4==tp4(2)) then
                                  arr2_tmp2(jj,7)=1;
								  arr2_tmp1(ii,4)=1;
							  endif
                          endif
                      endif
                  endif
				  if (tmp1==tp1(3)) then
                      if (tmp2==tp2(3)) then
                          if (tmp3==tp3(3)) then
                              if (tmp4==tp4(3)) then
                                  arr2_tmp2(jj,8)=1;
								  arr2_tmp1(ii,4)=1;
							  endif
                          endif
                      endif
                  endif
				  if (tmp1==tp1(4)) then
                      if (tmp2==tp2(4)) then
                          if (tmp3==tp3(4)) then
                              if (tmp4==tp4(4)) then
                                  arr2_tmp2(jj,9)=1;
								  arr2_tmp1(ii,4)=1;
							  endif
                          endif
                      endif
                  endif
			   enddo

			   enddo

               deallocate( arr1_tmp1,arr1_tmp2,arr1_tmp3,arr1_tmp4 );
			   deallocate( tp1,tp2,tp3,tp4 );


			   !============================================
			   !  Record the *.ipat
			   !============================================
		  tmp_filename2='.ipat'
                 call fstrcat (tmp_char,tmp_filename2,tmp_filename1)
                 open( Uinp,file=tmp_filename1 )
     
			   open(Uout_ipat,file=tmp_filename1)
			        write(Uout_ipat,*) maxfacepatch
					do jj=1,maxfacepatch
					   write(Uout_ipat,* ) arr2_tmp1(jj,1),arr2_tmp1(jj,2),&
					                          arr2_tmp1(jj,3),arr2_tmp1(jj,4)
					enddo
			   close(Uout_ipat);
               deallocate( arr2_tmp1 );

               write(*,*) "finish creating ",trim(tmp_filename1);

			   !============================================
			   !  Record the *.lv
			   !============================================
                                            tmp_filename2='.lv'
                 call fstrcat (tmp_char,tmp_filename2,tmp_filename1)
			   open(Uout_lv,file=tmp_filename1)
			        write(Uout_lv,*) maxlv
					do jj=1,maxlv
					   write(Uout_lv,*) arr2_tmp2(jj,1:9)
					enddo
			   close(Uout_lv);
               deallocate( arr2_tmp2 );

			   write(*,*) "finish creating ",trim(tmp_filename1);


			   
			    

		 
          
         


1201 format (3i8,6e16.9)
1202 format (3e20.7)
1203 format (19i8)
1204 format (3i16)
1205 format (9i12)
            end
