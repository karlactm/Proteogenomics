package Protein;

sub new {
    my $class = shift;
    my $self = {
        _id => shift,
        _description => shift,
        _sample => shift,
        _header => shift,   
        _sequence => shift
    };

    bless $self, $class;
    return $self;
}

sub getId {
    my $self = shift;
    return $self->{_id};
}

sub getDescription {
    my $self = shift;
    return $self->{_description};
}

sub getSample {
    my $self = shift;
    return $self->{_sample};
}

sub getHeader {
    my $self = shift;
    return $self->{_header};
}

sub getSequence {
    my $self = shift;
    return $self->{_sequence};
}

sub getSequenceSize {
    my $self = shift;
    return length $self->{_sequence};
}


1;
