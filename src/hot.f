* Written by Gilles Guillot february-august 2006
*
*
*


      subroutine hot(array,ngene,q,p,p1,p2,nset,res,
     &     m1,m2,v1,v2,w,bb,zz,ynum,index)
*     array : data
*     ngene : number of genes (number of rows)
*     q : number of genes used to discriminate (2,3 or 4)
*     p: number of arrays (number of columns)
*     p1,p2: nb of samples in each group (p1+p2=p)
*     nset: nb of sets of genes stored
*     res: matrix of results
*          nset lines , q+1 columns
*          q first columns give set the index of genes
*          last column give the T^2 statistics 
      implicit none 
      integer ngene,p,p1,p2,npair,q,nset
      double precision array,res

      integer g1,g2,g3,g4,g5,g6,j,j1,j2,index,indiv,iset,posleast,ii,
     &     INFO
      double precision m1,m2,v1,
     &     v2,w,mahal,f,least,bb,zz,
     &     YNUM,calcT2

      dimension array(ngene,p),res(nset,q+1),index(q),
     &     m1(q),m2(q),v1(q,q),v2(q,q),w(q,q),bb(q,q),zz(q),YNUM(q)

      write(*,*) 'debut de hot'
      
       
***********************
      if(q .eq. 2) then 
*     loop on pairs of genes
      iset = 0
      least = 1.d200
      do g1 = 1,ngene-1
         write(*,*) 'g1 =', g1
         do g2 = g1+1,ngene
            if(iset .le. nset) iset = iset + 1
c            write(*,*) ''
c            write(*,*) 'g1 =', g1
c            write(*,*) 'g2=',g2
c            write(*,*) 'iset=',iset
            index(1) = g1
            index(2) = g2
*     init of means and variances
            do j =1,q
               m1(j) = 0.
               m2(j) = 0.
            enddo
            do j1=1,q
               do j2=1,q
                  v1(j1,j2) = 0.  
                  v2(j1,j2) = 0.
                  w(j1,j2) = 0.
               enddo
            enddo
c$$$            write(*,*) 'm1=', m1
c$$$            write(*,*) 'm2=', m2
c$$$            write(*,*) 'v1=', v1
c$$$            write(*,*) 'v2=', v2
c$$$            write(*,*) ''

*     Means in group 1 and 2
            do j = 1,q
               do indiv = 1,p1
                  m1(j) = m1(j) + array(index(j),indiv)
               enddo
               do indiv = p1+1,p1+p2
                  m2(j) = m2(j) + array(index(j),indiv)
               enddo
               m1(j) = m1(j)/dble(p1)
               m2(j) = m2(j)/dble(p2)
            enddo

*     covariances in group 1 and 2
            do j1 = 1,q
               do j2 = j1,q
                  do indiv = 1,p1
                     v1(j1,j2) = v1(j1,j2) + 
     &                    (array(index(j1),indiv)-m1(j1))*
     &                    (array(index(j2),indiv)-m1(j2))
                  enddo
                  v1(j1,j2) = v1(j1,j2)/dble(p1-1)
                  do indiv = p1+1,p1+p2
                     v2(j1,j2) = v2(j1,j2) + 
     &                    (array(index(j1),indiv)-m2(j1))*
     &                    (array(index(j2),indiv)-m2(j2))
                  enddo
                  v2(j1,j2) = v2(j1,j2)/dble(p2-1)
                  w(j1,j2) = (p1*v1(j1,j2) + p2*v2(j1,j2))/(p1+p2)
                  w(j2,j1) = w(j1,j2) 
               enddo
            enddo

*     compute Fisher statistic (from Saporta p.351)
            f = ((dble(p1)+dble(p2)-dble(q)-1.)*dble(p1)*dble(p2)
     &           /(dble(q)*(dble(p1)+dble(p2))**2))*
     &           mahal(m1,m2,w,2,bb,YNUM,zz,INFO)
            if(INFO .ne. 0) then
               f = -999
            endif
 
c$$$            write(*,*) 'm1=', m1
c$$$            write(*,*) 'm2=', m2
c$$$            write(*,*) 'v1=', v1
c$$$            write(*,*) 'v2=', v2
c$$$            write(*,*) 'w=', w
c$$$c            write(*,*) 'mahal=', mahal(m1,m2,w,2)
c$$$            write(*,*) 'p1=',p1
c$$$            write(*,*) 'p2=',p2
c$$$            write(*,*) 'q=',q
c$$$            write(*,*) 'f=', f(iset)
c$$$            write(*,*) ''


*     store the result of the nset first sets
            if(iset .le. nset) then
               res(iset,1) = g1
               res(iset,2) = g2
               res(iset,3) = f
c               write(*,*) 'least=',least
               if(f .lt. least) then 
                  least = f
                  posleast = iset
c                  write(*,*) 'least=',least
c                  write(*,*) 'posleast=',posleast
               endif 
            endif
c$$$c            write(*,*) 'f=',f
c$$$c            write(*,*) 'least=',least
c$$$
*     now current set has to improve the criterion 
*     in order to be stored
            if(iset .gt. nset) then
               if(f .gt. least) then 
c                  write(*,*) 'coucou'
                  res(posleast,1) = g1
                  res(posleast,2) = g2
                  res(posleast,3) = f
c                  write(*,*) 'least avant=',least
c                  write(*,*) 'res(posleast,3)=',res(posleast,3)
                  least = 1.d200
                  do ii = 1,nset
                     if(res(ii,3) .lt. least) then 
                        least = res(ii,3)
                        posleast = ii
c                        write(*,*) 'res(ii,3)=',res(ii,3)
c                        write(*,*) 'least=',least
                     endif
                  enddo
c                  write(*,*) 'least apres=',least

               endif
            endif

         enddo
      enddo
      endif
*     endif (q .eq. 2)
***********************




***********************
      if(q .eq. 3) then 
      iset = 0
      least = 1.d200
      do g1 = 1,ngene-2
         write(*,*) 'g1 =', g1
         do g2 = g1+1,ngene-1
c             write(*,*) 'g2 =', g2
            do g3 = g2+1,ngene
c                write(*,*) 'g3 =', g3
               if(iset .le. nset) iset = iset + 1
c               write(*,*) 'iset=',iset
               index(1) = g1
               index(2) = g2
               index(3) = g3
               
*     init of means and variances
               do j =1,q
                  m1(j) = 0.
                  m2(j) = 0.
               enddo
               do j1=1,q
                  do j2=1,q
                     v1(j1,j2) = 0.
                     v2(j1,j2) = 0.
                     w(j1,j2) = 0.
                  enddo
               enddo
c                write(*,*) 'm1=', m1
c                write(*,*) 'm2=', m2
C                write(*,*) 'v1=', v1
C                write(*,*) 'v2=', v2
C                write(*,*) ''

*     Means in group 1 and 2
               do j = 1,q
                  do indiv = 1,p1
                     m1(j) = m1(j) + array(index(j),indiv)
                  enddo
                  do indiv = p1+1,p1+p2
                     m2(j) = m2(j) + array(index(j),indiv)
                  enddo
                  m1(j) = m1(j)/dble(p1)
                  m2(j) = m2(j)/dble(p2)
               enddo
c               write(*,*) 'm1=', m1
c               write(*,*) 'm2=', m2
*     covariances in group 1 and 2
               do j1 = 1,q
                  do j2 = j1,q
                     do indiv = 1,p1
                        v1(j1,j2) = v1(j1,j2) + 
     &                       (array(index(j1),indiv)-m1(j1))*
     &                       (array(index(j2),indiv)-m1(j2))
                     enddo
                     v1(j1,j2) = v1(j1,j2)/dble(p1-1)
                     do indiv = p1+1,p1+p2
                        v2(j1,j2) = v2(j1,j2) + 
     &                       (array(index(j1),indiv)-m2(j1))*
     &                       (array(index(j2),indiv)-m2(j2))
                     enddo
                     v2(j1,j2) = v2(j1,j2)/dble(p2-1)
                     w(j1,j2) = (p1*v1(j1,j2) + p2*v2(j1,j2))/
     &                    dble(p1+p2)
                     w(j2,j1) = w(j1,j2) 
                  enddo
               enddo
c               write(*,*) 'v1=', v1
c               write(*,*) 'v2=', v2
*     compute Fisher statistic (from Saporta p.351)
               f = ((dble(p1)+dble(p2)-dble(q)-1.)*
     &              dble(p1)*dble(p2)/
     &              (dble(q)*(dble(p1)+dble(p2))**2))*
     &              mahal(m1,m2,w,q,bb,YNUM,zz,INFO)
               if(INFO .ne. 0) then
                  f = -999
               endif
c$$$               write(*,*) 'g1 =', g1
c$$$               write(*,*) 'g2 =', g2
c$$$               write(*,*) 'g3 =', g3
c$$$               write(*,*) 'index=',index
c$$$               write(*,*) 'm1=', m1
c$$$               write(*,*) 'm2=', m2
c$$$               write(*,*) 'v1=', v1
c$$$               write(*,*) 'v2=', v2
c$$$               write(*,*) 'w=', w
c$$$               write(*,*) 'p1=',p1
c$$$               write(*,*) 'p2=',p2
c$$$               write(*,*) 'q=',q
c$$$               write(*,*) 'f=', f




c$$$               write(*,*) ''
c$$$
c$$$               write(*,*) 'iset=',iset
c$$$               write(*,*) 'nset=',nset
*     store the result of the nset first sets
               if(iset .le. nset) then
                  res(iset,1) = g1
                  res(iset,2) = g2
                  res(iset,3) = g3
                  res(iset,q+1) = f
c     write(*,*) 'least=',least
                  if(f .lt. least) then 
                     least = f
                     posleast = iset
c     write(*,*) 'least=',least
c                  write(*,*) 'posleast=',posleast
                  endif 
               endif


*     now current set has to improve the criterion 
*     in order to be stored
            if(iset .gt. nset) then
               if(f .gt. least) then 
c                  write(*,*) 'coucou'
                  res(posleast,1) = g1
                  res(posleast,2) = g2
                  res(posleast,3) = g3
                  res(posleast,q+1) = f
c                  write(*,*) 'least avant=',least
c                  write(*,*) 'res(posleast,q+1)=',res(posleast,q+1)
                  least = 1.d200
                  do ii = 1,nset
                     if(res(ii,q+1) .lt. least) then 
                        least = res(ii,q+1)
                        posleast = ii
c                        write(*,*) 'res(ii,q+1)=',res(ii,q+1)
c                        write(*,*) 'least=',least
                     endif
                  enddo
c                  write(*,*) 'least apres=',least

               endif
            endif

 

            enddo
c            write(*,*)  'end of loop on g3'
         enddo
c         write(*,*)  'end of loop on g2'
      enddo
c      write(*,*)  'end of loop on g1'
      endif
*     endif (q .eq. 3)
***********************


 



***********************
      if(q .eq. 4) then 
         iset = 0
         least = 1.d200
         do g1 = 1,ngene-3
c            write(*,*) 'g1 =', g1
            do g2 = g1+1,ngene-2
               write(*,*) 'g1,g2 =', g1,g2
               do g3 = g2+1,ngene-1
                  do g4 = g3+1,ngene
                     if(iset .le. nset) iset = iset + 1
                     index(1) = g1
                     index(2) = g2
                     index(3) = g3
                     index(4) = g4   

*     init of means and variances
                  do j =1,q
                     m1(j) = 0.
                     m2(j) = 0.
                  enddo
                  do j1=1,q
                     do j2=1,q
                        v1(j1,j2) = 0.
                        v2(j1,j2) = 0.
                        w(j1,j2) = 0.
                     enddo
                  enddo

*     Means in group 1 and 2
                     do j = 1,q
                        do indiv = 1,p1
                           m1(j) = m1(j) + array(index(j),indiv)
                        enddo
                        do indiv = p1+1,p1+p2
                           m2(j) = m2(j) + array(index(j),indiv)
                        enddo
                        m1(j) = m1(j)/dble(p1)
                        m2(j) = m2(j)/dble(p2)
                     enddo

*     covariances in group 1 and 2
                     do j1 = 1,q
                        do j2 = j1,q
                           do indiv = 1,p1
                              v1(j1,j2) = v1(j1,j2) + 
     &                             (array(index(j1),indiv)-m1(j1))*
     &                             (array(index(j2),indiv)-m1(j2))
                           enddo
                           v1(j1,j2) = v1(j1,j2)/dble(p1-1)
                           do indiv = p1+1,p1+p2
                              v2(j1,j2) = v2(j1,j2) + 
     &                             (array(index(j1),indiv)-m2(j1))*
     &                             (array(index(j2),indiv)-m2(j2))
                           enddo
                           v2(j1,j2) = v2(j1,j2)/dble(p2-1)
                           w(j1,j2) = (p1*v1(j1,j2) + p2*v2(j1,j2))/
     &                          dble(p1+p2)
                           w(j2,j1) = w(j1,j2) 
                        enddo
                     enddo

*     compute Fisher statistic (from Saporta p.351)
                     f = ((dble(p1)+dble(p2)-dble(q)-1.)*
     &                    dble(p1)*dble(p2)/
     &                    (dble(q)*(dble(p1)+dble(p2))**2))*
     &                    mahal(m1,m2,w,q,bb,YNUM,zz,INFO)
                     if(INFO .ne. 0) then
                        f = -999
                     endif

c$$$               write(*,*) 'g1 =', g1
c$$$               write(*,*) 'g2 =', g2
c$$$               write(*,*) 'g3 =', g3
c$$$               write(*,*) 'g4 =', g4
c$$$               write(*,*) 'index=',index
c$$$               write(*,*) 'm1=', m1
c$$$               write(*,*) 'm2=', m2
c$$$               write(*,*) 'v1=', v1
c$$$               write(*,*) 'v2=', v2
c$$$               write(*,*) 'w=', w
c$$$               write(*,*) 'p1=',p1
c$$$               write(*,*) 'p2=',p2
c$$$               write(*,*) 'q=',q
c$$$               write(*,*) 'f=', f
c$$$               write(*,*) ''
               
*     store the result of the nset first sets
                     if(iset .le. nset) then
                        res(iset,1) = g1
                        res(iset,2) = g2
                        res(iset,3) = g3
                        res(iset,4) = g4
                        res(iset,q+1) = f
                        if(f .lt. least) then 
                           least = f
                           posleast = iset
                        endif 
                     endif
                     
*     now current set has to improve the criterion 
*     in order to be stored
                     if(iset .gt. nset) then
                        if(f .gt. least) then 
                           res(posleast,1) = g1
                           res(posleast,2) = g2
                           res(posleast,3) = g3
                           res(posleast,4) = g4
                           res(posleast,q+1) = f
                           least = 1.d200
                           do ii = 1,nset
                              if(res(ii,q+1) .lt. least) then 
                                 least = res(ii,q+1)
                                 posleast = ii
                              endif
                           enddo
                        endif
                     endif
                     
                  enddo
c     write(*,*)  'end of loop on g4'
               enddo
c     write(*,*)  'end of loop on g3'
            enddo
c     write(*,*)  'end of loop on g2'
         enddo
c     write(*,*)  'end of loop on g1'
      endif
*     endif (q .eq. 4)
***********************





***********************
      if(q .eq. 5) then 
         iset = 0
         least = 1.d200
         do g1 = 1,ngene-4
c            write(*,*) 'g1 =', g1
            do g2 = g1+1,ngene-3
               write(*,*) 'g1,g2 =', g1,g2
               do g3 = g2+1,ngene-2
                  do g4 = g3+1,ngene-1
                     do g5 = g4+1,ngene
                        if(iset .le. nset) iset = iset + 1
                        index(1) = g1
                        index(2) = g2
                        index(3) = g3
                        index(4) = g4   
                        index(5) = g5
*     init of means and variances
                  do j =1,q
                     m1(j) = 0.
                     m2(j) = 0.
                  enddo
                  do j1=1,q
                     do j2=1,q
                        v1(j1,j2) = 0.
                        v2(j1,j2) = 0.
                        w(j1,j2) = 0.
                     enddo
                  enddo

*     Means in group 1 and 2
                     do j = 1,q
                        do indiv = 1,p1
                           m1(j) = m1(j) + array(index(j),indiv)
                        enddo
                        do indiv = p1+1,p1+p2
                           m2(j) = m2(j) + array(index(j),indiv)
                        enddo
                        m1(j) = m1(j)/dble(p1)
                        m2(j) = m2(j)/dble(p2)
                     enddo

*     covariances in group 1 and 2
                     do j1 = 1,q
                        do j2 = j1,q
                           do indiv = 1,p1
                              v1(j1,j2) = v1(j1,j2) + 
     &                             (array(index(j1),indiv)-m1(j1))*
     &                             (array(index(j2),indiv)-m1(j2))
                           enddo
                           v1(j1,j2) = v1(j1,j2)/dble(p1-1)
                           do indiv = p1+1,p1+p2
                              v2(j1,j2) = v2(j1,j2) + 
     &                             (array(index(j1),indiv)-m2(j1))*
     &                             (array(index(j2),indiv)-m2(j2))
                           enddo
                           v2(j1,j2) = v2(j1,j2)/dble(p2-1)
                           w(j1,j2) = (p1*v1(j1,j2) + p2*v2(j1,j2))/
     &                          dble(p1+p2)
                           w(j2,j1) = w(j1,j2) 
                        enddo
                     enddo

*     compute Fisher statistic (from Saporta p.351)
                     f = ((dble(p1)+dble(p2)-dble(q)-1.)*
     &                    dble(p1)*dble(p2)/
     &                    (dble(q)*(dble(p1)+dble(p2))**2))*
     &                    mahal(m1,m2,w,q,bb,YNUM,zz,INFO)
                     if(INFO .ne. 0) then
                        f = -999
                     endif

*     store the result of the nset first sets
                     if(iset .le. nset) then
                        res(iset,1) = g1
                        res(iset,2) = g2
                        res(iset,3) = g3
                        res(iset,4) = g4
                        res(iset,5) = g5
                        res(iset,q+1) = f
                        if(f .lt. least) then 
                           least = f
                           posleast = iset
                        endif 
                     endif

*     now current set has to improve the criterion 
*     in order to be stored
                     if(iset .gt. nset) then
                        if(f .gt. least) then 
                           res(posleast,1) = g1
                           res(posleast,2) = g2
                           res(posleast,3) = g3
                           res(posleast,4) = g4
                           res(posleast,5) = g5
                           res(posleast,q+1) = f
                           least = 1.d200
                           do ii = 1,nset
                              if(res(ii,q+1) .lt. least) then 
                                 least = res(ii,q+1)
                                 posleast = ii
                              endif
                           enddo
                        endif
                     endif
                     
                  enddo
c     write(*,*)  'end of loop on g5'
               enddo
c     write(*,*)  'end of loop on g4'
            enddo
c     write(*,*)  'end of loop on g3'
         enddo
c     write(*,*)  'end of loop on g2'
      enddo
c     write(*,*)  'end of loop on g1'
      endif                     
*     endif (q .eq. 5)
***********************


***********************
      if(q .eq. 6) then 
         iset = 0
         least = 1.d200
         do g1 = 1,ngene-5
c     write(*,*) 'g1 =', g1
            do g2 = g1+1,ngene-5
               write(*,*) 'g1,g2 =', g1,g2
               do g3 = g2+1,ngene-3
                  do g4 = g3+1,ngene-2
                     do g5 = g4+1,ngene-1
                        do g6 = g5+1,ngene
                           if(iset .le. nset) iset = iset + 1
                           index(1) = g1
                           index(2) = g2
                           index(3) = g3
                           index(4) = g4   
                           index(5) = g5
                           index(6) = g6
*     init of means and variances
                           do j =1,q
                              m1(j) = 0.
                              m2(j) = 0.
                           enddo
                           do j1=1,q
                              do j2=1,q
                                 v1(j1,j2) = 0.
                                 v2(j1,j2) = 0.
                                 w(j1,j2) = 0.
                              enddo
                           enddo
                           
*     Means in group 1 and 2
                           do j = 1,q
                              do indiv = 1,p1
                                 m1(j) = m1(j) + array(index(j),indiv)
                              enddo
                              do indiv = p1+1,p1+p2
                                 m2(j) = m2(j) + array(index(j),indiv)
                              enddo
                              m1(j) = m1(j)/dble(p1)
                              m2(j) = m2(j)/dble(p2)
                           enddo

*     covariances in group 1 and 2
                           do j1 = 1,q
                              do j2 = j1,q
                                 do indiv = 1,p1
                                    v1(j1,j2) = v1(j1,j2) + 
     &                                   (array(index(j1),indiv)-m1(j1))
     &                                  *(array(index(j2),indiv)-m1(j2))
                                 enddo
                                 v1(j1,j2) = v1(j1,j2)/dble(p1-1)
                                 do indiv = p1+1,p1+p2
                                    v2(j1,j2) = v2(j1,j2) + 
     &                                   (array(index(j1),indiv)-m2(j1))
     &                                  *(array(index(j2),indiv)-m2(j2))
                                 enddo
                                 v2(j1,j2) = v2(j1,j2)/dble(p2-1)
                                 w(j1,j2) = (p1*v1(j1,j2) + p2*v2(j1,j2)
     &                                )/dble(p1+p2)
                                 w(j2,j1) = w(j1,j2) 
                              enddo
                           enddo
                           
*     compute Fisher statistic (from Saporta p.351)
                           f = ((dble(p1)+dble(p2)-dble(q)-1.)*
     &                          dble(p1)*dble(p2)/
     &                          (dble(q)*(dble(p1)+dble(p2))**2))*
     &                          mahal(m1,m2,w,q,bb,YNUM,zz,INFO)
                           if(INFO .ne. 0) then
                              f = -999
                           endif
                           
*     store the result of the nset first sets
                           if(iset .le. nset) then
                              res(iset,1) = g1
                              res(iset,2) = g2
                              res(iset,3) = g3
                              res(iset,4) = g4
                              res(iset,5) = g5
                              res(iset,6) = g6
                              res(iset,q+1) = f
                              if(f .lt. least) then 
                                 least = f
                                 posleast = iset
                              endif 
                           endif
                           
*     now current set has to improve the criterion 
*     in order to be stored
                           if(iset .gt. nset) then
                              if(f .gt. least) then 
                                 res(posleast,1) = g1
                                 res(posleast,2) = g2
                                 res(posleast,3) = g3
                                 res(posleast,4) = g4
                                 res(posleast,5) = g5
                                 res(posleast,6) = g6
                                 res(posleast,q+1) = f
                                 least = 1.d200
                                 do ii = 1,nset
                                    if(res(ii,q+1) .lt. least) then 
                                       least = res(ii,q+1)
                                       posleast = ii
                                    endif
                                 enddo
                              endif
                           endif
                          
                        enddo
c     write(*,*)  'end of loop on g6'                        
                     enddo
c     write(*,*)  'end of loop on g5'
                  enddo
c     write(*,*)  'end of loop on g4'
               enddo
c     write(*,*)  'end of loop on g3'
            enddo
c     write(*,*)  'end of loop on g2'
         enddo
c     write(*,*)  'end of loop on g1'
      endif                     
*     endif (q .eq. 6)
***********************
      write(*,*) 'End of subroutine hot'
      end subroutine hot


 
  

***********************************************************
*
*     Mahalanobis distance (x-y)'A^-1(x-y)
*
      double precision function mahal(x,y,A,q,B,YNUM,z,INFO)
      implicit none
      integer q
      double precision x(q),y(q),A(q,q)
      integer i,j,INFO
      double precision B(q,q),d,z(q),YNUM(q)

c      write(*,*) 'x=',x
c      write(*,*) 'y=',y
      if(q .eq. 2) then
*     inverse of A
         d = A(1,1)*A(2,2) - A(1,2)*A(2,1)
         B(1,1) =  A(2,2)/d
         B(1,2) = -A(1,2)/d
         B(2,1) = -A(2,1)/d
         B(2,2) =  A(1,1)/d
         INFO = 0
      endif

      if(q .ge. 3) then
c         call solve3(A,B)
         call solve2(A,B,q,YNUM,INFO)
      endif

c      if(n .ge. 4) then
c         call solve(A,B,n,length)
c      endif
      
c      write(*,*) 'A=', A
c      write(*,*) 'B=', B
*     product z = A^-1(x-y)
      do i=1,q
         z(i) = 0.
         do j=1,q
            z(i) = z(i) + B(i,j)*(x(j)-y(j))
         enddo
      enddo
c      write(*,*) 'z=', z

*     product (x-y)*z
      mahal = 0.
      do i =1,q
         mahal = mahal + (x(i)-y(i))*z(i)
c         write(*,*) 'mahal=', mahal
      enddo
c      write(*,*) 'mahal=', mahal
      end
**************************************************************      



**************************************************************
*     invert a 3x3 matrix
      subroutine solve3(b,c)
      implicit none 
      double precision b(3,3), c(3,3)
      double precision det
      det  = (b(1,1)* b(2,2) - b(1,2)* b(2,1))* b(3,3) + 
     &     (b(1,3)* b(2,1) - b(1,1)* b(2,3))* b(3,2) + 
     &     (b(1,2)* b(2,3) - b(1,3)* b(2,2))* b(3,1)
      c(1,1) = (b(2,2)* b(3,3) - b(2,3)* b(3,2))/det
      c(2,1) = -(b(2,1)* b(3,3) - b(2,3)* b(3,1))/det
      c(3,1) = (b(2,1)* b(3,2) - b(2,2)* b(3,1))/det
	
      c(1,2) = -(b(1,2)* b(3,3) - b(1,3)* b(3,2))/det	
      c(2,2) = (b(1,1)* b(3,3) - b(1,3)* b(3,1))/det
      c(3,2) = -(b(1,1)* b(3,2) - b(1,2)* b(3,1))/det

      c(1,3) = (b(1,2)* b(2,3) - b(1,3)* b(2,2))/det
      c(2,3) = -(b(1,1)* b(2,3) - b(1,3)* b(2,1))/det
      c(3,3) = (b(1,1)* b(2,2) - b(1,2)* b(2,1))/det
      end subroutine solve3
************************************************************** 





**************************************************************
*     invert an nxn s.p.d matrix
      subroutine solve2(b,c,n,YNUM,INFO)
      implicit none 
      integer n
      double precision b(n,n), c(n,n)
      integer i,j,k


      INTEGER INFO,JOB
      DOUBLE PRECISION YNUM(n),RCOND,DET(2)
      parameter(JOB = 01)
      INFO = 0
*     cholesky decomposition 
      CALL dpoco (b, n, n, RCOND, Ynum, INFO)

*     inversion only ( JOB = 01 ): 
c      write(*,*) 'avant  dpodi'
      if(INFO .eq. 0) then
         CALL dpodi (b, n, n ,DET, JOB)
         do j=1,n
            do i=1,j
               c(i,j) = b(i,j)
               c(j,i) = c(i,j)
            enddo
         enddo
      endif
         
c      write(*,*) 'fin de solve2'
      end subroutine solve2
********************************************************


