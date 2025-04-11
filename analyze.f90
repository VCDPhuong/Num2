program extract_lammpstrj
    implicit none
    character(len=256) :: line
    integer :: i, num_atoms, atom_id, timestep
    real :: x, y, z, vx, vy, vz
    logical :: is_position, is_velocity
    character(len=256) :: filename
    character(len=256) :: output_file
    integer :: input_unit, output_unit

    ! Input and output file names
    filename = "trajectory.lammpstrj"
    output_file = "extracted_data.txt"

    ! Open the input and output files
    open(newunit=input_unit, file=filename, status='old', action='read')
    open(newunit=output_unit, file=output_file, status='replace', action='write')

    is_position = .false.
    is_velocity = .false.
    timestep = -1

    do
        read(input_unit, '(A)', iostat=i) line
        if (i /= 0) exit  ! Exit loop at end of file

        ! Check for timestep
        if (index(line, "ITEM: TIMESTEP") > 0) then
            read(input_unit, '(A)', iostat=i) line
            read(line, *) timestep
            write(output_unit, '(A, I8)') "Timestep: ", timestep
            cycle
        end if

        ! Check for position or velocity section
        if (index(line, "ITEM: ATOMS") > 0) then
            if (index(line, "x y z") > 0) then
                is_position = .true.
                is_velocity = .false.
                write(output_unit, '(A)') "Position Data:"
            else if (index(line, "vx vy vz") > 0) then
                is_position = .false.
                is_velocity = .true.
                write(output_unit, '(A)') "Velocity Data:"
            else
                is_position = .false.
                is_velocity = .false.
            end if
            cycle
        end if

        ! Extract data if in position or velocity section
        if (is_position .or. is_velocity) then
            read(line, *) atom_id, x, y, z, vx, vy, vz
            if (is_position) then
                write(output_unit, '(I8, 3F12.6)') atom_id, x, y, z
            else if (is_velocity) then
                write(output_unit, '(I8, 3F12.6)') atom_id, vx, vy, vz
            end if
        end if
    end do

    ! Close files
    close(input_unit)
    close(output_unit)

    print *, "Data extraction complete. Output written to ", output_file
end program extract_lammpstrj