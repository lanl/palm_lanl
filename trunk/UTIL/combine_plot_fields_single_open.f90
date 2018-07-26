 PROGRAM combine_plot_fields

!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2018  Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! ---------------------
! $Id: combine_plot_fields_single_open.f90 2718 2018-01-02 08:49:38Z maronga $
! Corrected "Former revisions" section
! 
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
!
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores (not tested)
!
! code put under GPL (PALM 3.9)
! Prozessordateien werden einzeln geoeffnet und geschlossen
!
! 23/02/99  Keine Bearbeitung komprimierter 3D-Daten
! Ursprungsversion vom 28/07/97
!
!
! Description:
! -------------
! Vereinigung der von PARLES im Parallelbetrieb erzeugten Teilfelder zu
! gemeinsamen, das jeweilige Gesamtgebiet umfassenden Plotfeldern
!------------------------------------------------------------------------------!

    IMPLICIT NONE

!
!-- Lokale Variablen
    CHARACTER (LEN=2) ::  modus
    CHARACTER (LEN=9) ::  id_char

    INTEGER, PARAMETER ::  spk = SELECTED_REAL_KIND( 6 )

    INTEGER ::  danz, fanz, id, nxa, nxag, nxe, nxeg, nya, nyag, nye, nyeg, &
                nza, nzag, nze, nzeg, xa, xe, ya, ye, za, ze

    LOGICAL ::  compressed, found

    REAL ::  dx
    REAL, DIMENSION(:),   ALLOCATABLE   ::  eta, ho, hu
    REAL, DIMENSION(:,:), ALLOCATABLE   ::  pf
    REAL(spk), DIMENSION(:,:,:), ALLOCATABLE ::  pf3d


!
!-- 2D-Felder fuer ISO2D
!-- Hauptschleife ueber die 3 Schnittarten, beginnend mit xy-Schnitt
    modus = 'XY'
    PRINT*, ''
    DO  WHILE ( modus == 'XY'  .OR.  modus == 'XZ'  .OR.  modus == 'YZ' )
!
!--    Pruefen, ob Basisdatei von PE0 vorhanden
       danz = 0
       WRITE (id_char,'(A2,''_'',I6.6)')  modus, danz
       INQUIRE ( FILE='PLOT2D_'//id_char, EXIST=found )
!
!--    Anzahl der Dateien feststellen
       DO  WHILE ( found )

          danz = danz + 1
          WRITE (id_char,'(A2,''_'',I6.6)')  modus, danz
          INQUIRE ( FILE='PLOT2D_'//id_char, EXIST=found )

       ENDDO
!
!--    Info-Ausgabe
       PRINT*, ''
       PRINT*, '*** combine_plot_fields ***'
       IF ( danz /= 0 )  THEN
          PRINT*, modus,'-Schnitt:  ', danz, ' Datei(en) gefunden'
       ELSE
          PRINT*, 'keine ', modus, '-Schnitte vorhanden'
       ENDIF

!
!--    Einzelne Felder einlesen, bis keine mehr auf den Dateien vorhanden
       fanz = 0
       DO  WHILE ( danz /= 0 )

!
!--       Schleife ueber alle Dateien
          DO  id = 0, danz-1
!
!--          Prozessordatei oeffnen
             WRITE (id_char,'(A2,''_'',I6.6)')  modus, id
             OPEN ( 1, FILE='PLOT2D_'//id_char, FORM='UNFORMATTED', &
                       POSITION='ASIS' )

!
!--          Bei erster Datei und erstem eingelesenen Teilfeld Gesamtfeld
!--          allokieren, Ausgabedatei oeffnen und Koordinateninformationen
!--          schreiben
             IF ( id == 0  .AND.  fanz == 0 )  THEN
                READ ( 1 )  nxag, nxeg, nyag, nyeg
                READ ( 1 )  nxa, nxe, nya, nye
                ALLOCATE ( eta(nya:nye), ho(nxa:nxe), hu(nxa:nxe), &
                           pf(nxag:nxeg,nyag:nyeg) )
                READ ( 1 )  dx, eta, hu, ho

                OPEN ( 2, FILE='PLOT2D_'//modus, FORM='UNFORMATTED' )
                WRITE ( 2 )  dx, eta, hu, ho
             ENDIF
!
!--          Teilfeld einlesen und ausgeben
             READ ( 1, END=998 )  xa, xe, ya, ye
!
!--          Falls PE kein Teilfeld ausgegeben hat, sind Indices entsprechend
!--          gesetzt
             IF ( .NOT. ( xa == -1  .AND.  xe == -1  .AND. &
                          ya == -1  .AND.  ye == -1 ) )  THEN
                READ ( 1 )  pf(xa:xe,ya:ye)
             ENDIF
             IF ( id == 0 )  fanz = fanz + 1

!
!--          Prozessordatei schliessen
             CLOSE ( 1 )

          ENDDO
!
!--       Ausgabe des jeweiligen Gesamtfeldes
          WRITE ( 2 )  pf(nxa:nxe,nya:nye)
       
       ENDDO

998    IF ( danz /= 0 )  THEN
!
!--       Anzahl der bearbeiteten Felder ausgeben
          PRINT*, modus, '-Schnitt:  ', fanz, ' Feld(er) ausgegeben'
!
!--       Dateien schliessen, allokierte Felder freigeben
          CLOSE ( 1 )
          CLOSE ( 2 )
          DEALLOCATE ( eta, ho, hu, pf )
       ENDIF
!
!--    Naechste Schnittebene
       SELECT CASE ( modus )
          CASE ( 'XY' )
             modus = 'XZ'
          CASE ( 'XZ' )
             modus = 'YZ'
          CASE ( 'YZ' )
             modus = 'no'
       END SELECT

    ENDDO


!
!-- 3D-Felder fuer AVS
!
!-- Pruefen, ob Basisdatei von PE0 vorhanden
    danz = 0
    WRITE (id_char,'(I6.6)')  danz
    INQUIRE ( FILE='PLOT3D_DATA_'//TRIM( id_char ), EXIST=found )

!
!-- Vereinigung darf nur erfolgen, wenn 3D-Daten unkomprimiert vorliegen
    INQUIRE ( FILE='PLOT3D_COMPRESSED', EXIST=compressed )

!
!-- Anzahl der Dateien feststellen
    DO  WHILE ( found  .AND.  .NOT. compressed )

       danz = danz + 1
       WRITE (id_char,'(I6.6)')  danz
       INQUIRE ( FILE='PLOT3D_DATA_'//TRIM( id_char ), EXIST=found )

    ENDDO

!
!-- Info-Ausgabe
    PRINT*, ' '
    PRINT*, '*** combine_plot_fields ***'
    IF ( danz /= 0 )  THEN
       PRINT*, '3D-Ausgabe:  ', danz, ' Datei(en) gefunden'
    ELSE
       IF ( found .AND. compressed )  THEN
          PRINT*, '3D-Ausgabe nicht vorgenommen, da Daten komprimiert vorliegen'
       ELSE
          PRINT*, 'keine 3D-Ausgaben vorhanden'
       ENDIF
    ENDIF

!
!-- Einzelne Felder einlesen, bis keine mehr auf den Dateien vorhanden
    fanz = 0
    DO  WHILE ( danz /= 0 )

!
!--    Schleife ueber alle Dateien
       DO  id = 0, danz-1
!
!--       Prozessordatei oeffnen
          WRITE (id_char,'(I6.6)')  id
          OPEN ( 1, FILE='PLOT3D_DATA_'//TRIM( id_char ), FORM='UNFORMATTED', &
                    POSITION='ASIS' )
!
!--       Bei erster Datei und erstem eingelesenen Teilfeld Gesamtfeld
!--       allokieren und Ausgabedatei oeffnen
          IF ( id == 0  .AND.  fanz == 0 )  THEN
             READ ( 1 )  nxag, nxeg, nyag, nyeg, nzag, nzeg
             READ ( 1 )  nxa, nxe, nya, nye, nza, nze
             ALLOCATE ( pf3d(nxag:nxeg,nyag:nyeg,nzag:nzeg) )
             OPEN ( 2, FILE='PLOT3D_DATA', FORM='UNFORMATTED' )
          ENDIF
!
!--       Teilfeld einlesen und ausgeben
          READ ( 1, END=999 )  xa, xe, ya, ye, za, ze
          READ ( 1 )  pf3d(xa:xe,ya:ye,za:ze)
          IF ( id == 0 )  fanz = fanz + 1

!
!--       Prozessordatei schliessen
          CLOSE ( 1 )

       ENDDO
!
!--    Ausgabe des jeweiligen Gesamtfeldes
       WRITE ( 2 )  pf3d(nxa:nxe,nya:nye,nza:nze)
       
    ENDDO

999 IF ( danz /= 0 )  THEN
!
!--    Anzahl der bearbeiteten Felder ausgeben
       PRINT*, '3D-Ausgabe:  ', fanz, ' Feld(er) ausgegeben'
!
!--    Dateien schliessen, allokierte Felder freigeben
       CLOSE ( 1 )
       CLOSE ( 2 )
       DEALLOCATE ( pf3d )
    ENDIF

 END PROGRAM combine_plot_fields
