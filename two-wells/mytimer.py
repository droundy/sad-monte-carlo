import time

class Timer: 
    # Initializing
    def __init__(self, message):
        self.start = time.process_time()
        self.message = message
  
    # Calling destructor
    def __del__(self):
        t = time.process_time() - self.start
        if t > 60:
            print(self.message, 'took %.2g minutes' % (t/60))
        else:
            print(self.message, 'took %.2g sec' % t)
