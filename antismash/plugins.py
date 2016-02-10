import straight.plugin


def load_plugins():
    '''Load all available types of plugins'''
    plugins = {
        'cluster': [],
        'generic': [],
        'genome_wide': [],
        'input': [],
        'output': [],
        'specific': []
    }

    for key in plugins.keys():
        plugins[key] = list(straight.plugin.load('antismash.{}'.format(key)))
        plugins[key].sort(key=lambda x: x.priority)

    return plugins


def check_prereqs(plugins, options):
    '''Check if all prerequisites are met'''
    return []
