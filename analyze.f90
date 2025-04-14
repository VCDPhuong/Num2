
program main
 implicit none
    integer :: num_atoms
    real(8), allocatable :: positions(:,:), velocities(:,:), timestep(:)
    character(len=256) :: filename, output_filename
    integer :: i, j
    allocate(positions(3, num_atoms))
    allocate(velocities(3, num_atoms))
    ! Open the file
    open(10, file='Data/Output.lammpstrj', action='read', iostat=i)
    if (i /= 0) then
        print*, "Error: Unable to open file 'Data/Output.lammpstrj'"
        stop
    end if
    read(10, '(A)', iostat=i) filename
    read
    print*, "File name: ", filename
    !call load_velocity_position(num_atoms, positions, velocities)
    close(10)
    print*, positions(:,1), velocities(:,1)
    contains
subroutine load_velocity_position(num_atoms, positions, velocities)
    implicit none
    integer, intent(out) :: num_atoms
    real(8), allocatable, intent(out) :: positions(:,:), velocities(:,:)
    integer :: i
    character(len=256) :: line
    real(8) :: x, y, z, vx, vy, vz
    integer :: id


    ! Read the file line by line
    do
        read(10, '(A)', iostat=i) line
        if (i /= 0) exit

        ! Check for specific ITEM lines
        if (trim(line) == 'ITEM: NUMBER OF ATOMS') then
            read(10, *) num_atoms
            print *, num_atoms
            allocate(positions(3, num_atoms))
            allocate(velocities(3, num_atoms))
        else if (trim(line) == 'ITEM: ATOMS id x y z vx vy vz') then
            do i = 1, num_atoms
                read(10, *) id, x, y, z, vx, vy, vz
                positions(:, id) = (/ x, y, z /)
                velocities(:, id) = (/ vx, vy, vz /)
            end do
        end if
    end do

    ! Close the file
end subroutine load_velocity_position
end program