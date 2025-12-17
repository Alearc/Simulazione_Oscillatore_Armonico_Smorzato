!************************Modulo contente i kind utilizzati nel programma*********************************************
module prec 
implicit none
   integer,parameter :: rk=selected_real_kind(16)  !kind per i real
   integer, parameter :: ik=selected_int_kind(12)  !kind per gli integer
end module prec
!*********************************************************************************************************************

!************************INIZIO PROGRAMMA*****************************************************************************
program harm_smorz
 use prec                              !il programma utilizza il modulo prec
 implicit none
 real(kind=rk)      :: massa,kappa,c   ! valori di default per massa e cost. elastica
 real(kind=rk)      :: dt,ekin,epot,pos0 
 real(kind=rk)      :: ekin_an,epot_an,gam,w
 real(kind=rk)      :: pos,vel,vel_parziale,f,fp,ff
 real(kind=rk)      :: pos_an,vel_an,pos0_an
 integer(kind=ik)   :: nstep,it
 write(unit=*,fmt="(a)",advance="no")"delta t : "    
 read*,dt                                            
 write(unit=*,fmt="(a)",advance="no")"n.step: "      
 read*,nstep                                        
 write(unit=*,fmt="(a)",advance="no")"massa: "     
 read*,massa
 write(unit=*,fmt="(a)",advance="no")"kappa: "
 read*,kappa
 write(unit=*,fmt="(a)",advance="no")"c, costante di attrito lineare:"
 read*,c
 write(unit=*,fmt="(a)",advance="no")"pos(0): "
 read*,pos
 write(unit=*,fmt="(a)",advance="no")"vel(0): "
 read*,vel
 write(unit=*,fmt="(a)",advance="no")"pos(0), formula analitica: "
 read*,pos_an
 write(unit=*,fmt="(a)",advance="no")"vel(0), formula analitica: "
 read*,vel_an
 

 it=0_ik                                          ! step 0 : valori iniziali 
 write(unit=1,fmt=*)it,it*dt,pos,vel
 epot =  0.5 * kappa * pos**2                     !energia potenziale elastica
 f    = - kappa * pos - c * vel                   !legge di Hooke e legge di Stokes
 ekin =  0.5 * massa * vel**2                     !energia cinetica della massa
 write(unit=2,fmt=*)it,dt*it,ekin,epot,ekin+epot  !scrittura su file dei valori iniziali che verranno poi sostituiti
                                                  !da quelli ottenuti con l'algoritmo

 !***********************parte relativa alla formula analitica****************************************************************
 gam = c/(2*massa)                                !costante per la formula analitica
 w = sqrt(kappa/massa - (c/(2*massa))**2)         !pulsazione
 pos_an = pos_an*exp(-gam*it*dt)*cos(w*it*dt)+(vel_an+gam*pos_an)/w*exp(-gam*it*dt)*sin(w*it*dt)         ! posizione analitica
 vel_an = vel_an*cos(w*it*dt)*exp(-gam*it*dt)-(pos_an*gam**2 +gam*vel_an + pos_an*w**2)/w * sin(w*it*dt) !velocità analitica
 
 ekin_an = 0.5 * massa * vel_an**2                !energia cinetica analitica
 epot_an =  0.5 * kappa * pos_an**2               !energia potenziale analitica
 write(unit=4,fmt=*)it*dt,pos_an,vel_an           !scrittura su file dei valori iniziali che verranno poi sostituiti 
 write(unit=3,fmt=*)it*dt,"ekin:",ekin_an,"epot:",epot_an,"en_tot:", ekin_an+epot_an  !da quelli ottenuti con la formula
 !*****************************************************************************************************************************


 !*************INIZIO IMPLEMENTAZIONE ALGORITMO EULER-TRAPEZOIDALE***************************************************************
 do it = 1,nstep
    pos0 = pos                                              !variabile temporanea necessaria per la formula iniziale della forza
    f = - kappa * pos0 - c * vel                            !che contiene sempre il valore iniziale della posizione in ogni ciclo.
    pos = pos + vel * dt + 0.5 * f/massa * dt**2            !valore finale della posizione con algoritmo
    vel_parziale = vel + f/massa * dt                       !valore parziale della velocità con algoritmo
    fp = - kappa * pos - c * vel_parziale                   !valore paziale della forza con algotitmo
    epot =  0.5 * kappa * pos**2                            !energia potenziale con algoritmo
    vel = vel + 0.5 * (f/massa + fp/massa) * dt             !valore finale della velocità con algoritmo
    ff = - kappa * pos - c * vel                            !valore finale della forza con algoritmo
    write(unit=1,fmt=*)it,it*dt,pos,vel
    ekin = 0.5 * massa * vel**2                             !energia cinetica con algoritmo
    write(unit=2,fmt=*)it,it*dt,"ekin:",ekin,"epot:",epot,"en_tot:", ekin+epot
 !********************************************************************************************************************************** 
    
 !*************************parte relativa alla formula analitica********************************************************************
    pos0_an = pos_an                                        !variabile temporanea necessaria per la formula analitica della velocità
    pos_an = pos_an*exp(-gam*dt)*cos(w*dt)+(vel_an+gam*pos_an)/w*exp(-gam*dt)*sin(w*dt)  !posizione finale con formula analitica
    vel_an = vel_an*cos(w*dt)*exp(-gam*dt)-(pos0_an*gam**2 +gam*vel_an + pos0_an*w**2)/w * sin(w*dt) !velocità fin. con form analitica
    write(unit=4,fmt=*)dt*it,pos_an,vel_an
    ekin_an = 0.5 * massa * vel_an**2                       !energia cinetica con formula analitica
    epot_an =  0.5 * kappa * pos_an**2                      !energia potenziale con formula analitica
    write(unit=3,fmt=*)dt*it,"ekin_an:",ekin_an,"epot_an:",epot_an,"en_tot_an:", ekin_an+epot_an
 !**********************************************************************************************************************************   
 end do
end program harm_smorz

