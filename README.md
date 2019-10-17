# mr.ash

The mr.ash R package implements mr.ash

## Quick Start

To install the latest version of varbvs2 from GitHub, clone or
download the git repository, then use the `install_local` function
from [devtools][devtools]. Assuming your working directory contains
the varbvs2 repository, run this code to install the package:

```R
install.packages("devtools")
library(devtools)
list.files(pattern = "mr.ash") # Should output "varbvs".
install_local("mr.ash")
```

This command should automatically install all required packages if
they are not installed already.


## License

Copyright (c) 2018-2019, Youngseok Kim, Peter Carbonetto and Matthew
Stephens.

All source code and software in this repository are made available
under the terms of the [MIT license][mit-license]. See
file [LICENSE](LICENSE) for the full text of the license.

## Credits

The mr.ash R package was developed by [Youngseok Kim][youngseok] and
[Peter Carbonetto][peter] at the [University of Chicago][uchicago],
with help from [Matthew Stephens][matthew].

[mit-license]: https://opensource.org/licenses/mit-license.html
[devtools]: https://github.com/r-lib/devtools
[uchicago]: https://www.uchicago.edu
[youngseok]: https://github.com/youngseok-kim
[peter]: https://pcarbo.github.io
[matthew]: http://stephenslab.uchicago.edu
