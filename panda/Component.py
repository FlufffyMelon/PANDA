class Component:
    def __init__(
        self, name, region, insertion_region=None, density=None, numbers=None, **kwargs
    ):
        self.name = name
        self.region = region
        self.insertion_region = (
            insertion_region if insertion_region is not None else region
        )

        region_volume = self.region.get_volume()
        if density is not None and numbers is None:
            self.density = density
            self.numbers = int(round(region_volume * density))
        elif numbers is not None and density is None:
            self.numbers = numbers
            self.density = numbers / region_volume if region_volume > 0 else 0
        elif density is not None and numbers is not None:
            self.density = density
            self.numbers = numbers
        else:
            raise ValueError(
                "Either density or numbers must be specified for Component {}".format(
                    name
                )
            )
